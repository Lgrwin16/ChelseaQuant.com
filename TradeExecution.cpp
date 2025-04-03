#include "TradeExecution.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <numeric>

TradeExecution::TradeExecution(
    std::shared_ptr<IBConnector> ibConnector,
    std::shared_ptr<Statistics> statistics,
    double initialCapital,
    const PositionSizing& positionSizing)
    : ibConnector(ibConnector),
      statistics(statistics),
      positionSizing(positionSizing),
      accountCapital(initialCapital),
      availableCapital(initialCapital),
      totalPortfolioValue(initialCapital),
      totalRiskExposure(0.0) {

    std::cout << "Trade execution module initialized with $"
              << initialCapital << " of initial capital" << std::endl;
}

TradeExecution::~TradeExecution() {
    // Clean up resources if needed
}

void TradeExecution::updateAccountInfo() {
    if (!ibConnector->isConnected()) {
        std::cerr << "Error: Not connected to Interactive Brokers" << std::endl;
        return;
    }

    // In a real implementation, this would fetch account info from IB API
    accountCapital = ibConnector->getAccountCash();

    // Update position information
    currentPositions.clear();
    positionValues.clear();

    std::vector<std::string> positions = ibConnector->getPositions();
    for (const auto& symbol : positions) {
        double posSize = ibConnector->getPositionSize(symbol);
        currentPositions[symbol] = static_cast<int>(posSize);

        double currentPrice = ibConnector->getLatestPrice(symbol).close;
        positionValues[symbol] = std::abs(posSize) * currentPrice;
    }

    // Calculate total portfolio value and available capital
    double positionsTotal = std::accumulate(
        positionValues.begin(), positionValues.end(), 0.0,
        [](double sum, const auto& pair) { return sum + pair.second; }
    );

    totalPortfolioValue = accountCapital + positionsTotal;
    availableCapital = accountCapital;

    // Calculate total risk exposure
    updatePortfolioMetrics();

    std::cout << "Account info updated: Total portfolio value: $"
              << totalPortfolioValue << ", Available capital: $"
              << availableCapital << std::endl;
}

void TradeExecution::updatePortfolioMetrics() {
    totalRiskExposure = 0.0;

    for (const auto& [symbol, quantity] : currentPositions) {
        if (quantity == 0) continue;

        double positionValue = positionValues[symbol];

        // Calculate position risk based on volatility
        PriceTimeSeries priceData = ibConnector->getHistoricalData(
            symbol, "TODAY", 30, "1 day", "ADJUSTED_LAST"
        );

        std::vector<double> prices;
        for (const auto& data : priceData) {
            prices.push_back(data.close);
        }

        double volatility = statistics->calculateVolatility(prices);
        double beta = 1.0;  // Default beta for single asset

        try {
            beta = statistics->calculateBeta(symbol);
        } catch (const std::exception& e) {
            // If beta calculation fails, use default
        }

        // VaR-based risk contribution (simplified)
        double positionRisk = positionValue * volatility * std::sqrt(10) / 100.0;  // 10-day VaR at ~2 stddev
        totalRiskExposure += positionRisk;
    }

    std::cout << "Portfolio heat: " << calculatePortfolioHeat() * 100.0 << "%" << std::endl;
}

double TradeExecution::getAccountCapital() const {
    return accountCapital;
}

double TradeExecution::getAvailableCapital() const {
    return availableCapital;
}

double TradeExecution::getTotalPortfolioValue() const {
    return totalPortfolioValue;
}

double TradeExecution::getTotalRiskExposure() const {
    return totalRiskExposure;
}

double TradeExecution::calculatePortfolioHeat() const {
    if (totalPortfolioValue <= 0.0) {
        return 0.0;
    }
    return totalRiskExposure / totalPortfolioValue;
}

int TradeExecution::calculateShares(const TradeSignal& signal) {
    double currentPrice = ibConnector->getLatestPrice(signal.symbol).close;
    if (currentPrice <= 0.0) {
        std::cerr << "Error: Invalid price for " << signal.symbol << std::endl;
        return 0;
    }

    double positionSizeUSD = 0.0;

    switch (positionSizing.model) {
        case PositionSizingModel::FIXED_PERCENTAGE:
            // Simple fixed percentage of total capital
            positionSizeUSD = totalPortfolioValue * positionSizing.maxPositionSize;
            break;

        case PositionSizingModel::VOLATILITY_ADJUSTED:
            positionSizeUSD = calculateVolatilityAdjustedSize(signal);
            break;

        case PositionSizingModel::KELLY_CRITERION:
            positionSizeUSD = calculateKellyCriterion(signal) * totalPortfolioValue;
            // Usually we use half-Kelly for safety
            positionSizeUSD *= 0.5;
            break;

        case PositionSizingModel::EQUAL_RISK:
            // Allocate equal risk per position
            {
                // Base position size on maximum risk per trade
                double riskAmount = totalPortfolioValue * positionSizing.riskPerTrade;

                // Calculate stop loss distance in percentage
                double stopLossPct = positionSizing.stopLossPct;
                if (signal.stopLossPrice > 0) {
                    stopLossPct = std::abs(currentPrice - signal.stopLossPrice) / currentPrice;
                }

                // Size position such that loss at stop = riskAmount
                positionSizeUSD = riskAmount / stopLossPct;
            }
            break;

        case PositionSizingModel::VALUE_AT_RISK:
            // VaR-based position sizing
            {
                double varLimit = totalPortfolioValue * positionSizing.riskPerTrade;
                double var = calculateValueAtRisk(signal, 1, 0.95); // 1-day 95% VaR per dollar

                positionSizeUSD = varLimit / var;
            }
            break;
    }

    // Apply maximum position size constraint
    double maxPositionSizeUSD = totalPortfolioValue * positionSizing.maxPositionSize;
    positionSizeUSD = std::min(positionSizeUSD, maxPositionSizeUSD);

    // Apply available capital constraint
    positionSizeUSD = std::min(positionSizeUSD, availableCapital);

    // Calculate number of shares
    int shares = static_cast<int>(positionSizeUSD / currentPrice);

    // Ensure minimum trade size
    if (shares < 1) {
        shares = 0;
    }

    // Log position size calculation
    std::cout << "Position sizing for " << signal.symbol
              << ": $" << positionSizeUSD
              << " (" << shares << " shares at $" << currentPrice << ")" << std::endl;

    return shares;
}

std::pair<int, int> TradeExecution::calculatePairShares(
    const TradeSignal& signal, const TradeSignal& hedgeSignal) {

    // Get current prices
    double priceMain = ibConnector->getLatestPrice(signal.symbol).close;
    double priceHedge = ibConnector->getLatestPrice(hedgeSignal.symbol).close;

    if (priceMain <= 0.0 || priceHedge <= 0.0) {
        std::cerr << "Error: Invalid prices for pair" << std::endl;
        return {0, 0};
    }

    // Calculate hedge ratio if not provided
    double hedgeRatio = signal.hedgeRatio;
    if (hedgeRatio <= 0.0) {
        hedgeRatio = statistics->optimalHedgeRatio(signal.symbol, hedgeSignal.symbol);
    }

    // Calculate total position size based on portfolio constraints
    double totalPairExposure = totalPortfolioValue * positionSizing.maxPositionSize;

    // Split exposure between the pair
    double mainExposure = totalPairExposure / (1.0 + hedgeRatio * priceHedge / priceMain);
    double hedgeExposure = totalPairExposure - mainExposure;

    // Apply available capital constraint
    double availableExposure = std::min(totalPairExposure, availableCapital);
    mainExposure = mainExposure * (availableExposure / totalPairExposure);
    hedgeExposure = hedgeExposure * (availableExposure / totalPairExposure);

    // Calculate number of shares
    int sharesMain = static_cast<int>(mainExposure / priceMain);
    int sharesHedge = static_cast<int>(hedgeExposure / priceHedge);

    // Adjust to maintain proper hedge ratio
    if (sharesMain > 0 && sharesHedge > 0) {
        // Calculate actual hedge ratio with integer shares
        double actualRatio = (sharesHedge * priceHedge) / (sharesMain * priceMain);

        // If actual ratio deviates too much, adjust shares
        if (std::abs(actualRatio - hedgeRatio) / hedgeRatio > 0.05) {
            // Adjust hedge shares to match the target ratio more closely
            sharesHedge = static_cast<int>(sharesMain * priceMain * hedgeRatio / priceHedge);
        }
    }

    std::cout << "Pair position sizing: " << signal.symbol << " - " << hedgeSignal.symbol
              << " with hedge ratio " << hedgeRatio << std::endl;
    std::cout << "  " << signal.symbol << ": " << sharesMain << " shares ($" << mainExposure << ")" << std::endl;
    std::cout << "  " << hedgeSignal.symbol << ": " << sharesHedge << " shares ($" << hedgeExposure << ")" << std::endl;

    return {sharesMain, sharesHedge};
}

double TradeExecution::calculateVolatilityAdjustedSize(const TradeSignal& signal) const {
    // Calculate position size inversely proportional to volatility
    double volatility = signal.volatility;
    if (volatility <= 0.0) {
        // If volatility not provided, use historical
        PriceTimeSeries priceData = ibConnector->getHistoricalData(
            signal.symbol, "TODAY", 30, "1 day", "ADJUSTED_LAST"
        );

        std::vector<double> prices;
        for (const auto& data : priceData) {
            prices.push_back(data.close);
        }

        volatility = statistics->calculateVolatility(prices);
    }

    if (volatility <= 0.001) {
        volatility = 0.001;  // Minimum volatility to avoid division by zero
    }

    // Base position size on risk per trade
    double riskAmount = totalPortfolioValue * positionSizing.riskPerTrade;

    // Adjust position size based on volatility (inverse relationship)
    // Target volatility of 20% annualized as reference
    double targetVol = 0.20;
    double scalingFactor = targetVol / volatility;

    // Apply reasonable bounds to the scaling factor
    scalingFactor = std::max(0.2, std::min(scalingFactor, 5.0));

    double positionSizeUSD = riskAmount * scalingFactor * 10.0;  // 10x leverage on risk amount

    return positionSizeUSD;
}

double TradeExecution::calculateKellyCriterion(const TradeSignal& signal) const {
    // Kelly Criterion: f* = (p * b - (1 - p)) / b
    // where p = probability of winning, b = win/loss ratio

    // For simplicity, we'll use confidence score as a proxy for p
    double winProbability = 0.5 + (signal.confidenceScore / 10.0);  // Convert confidence to probability
    winProbability = std::max(0.1, std::min(winProbability, 0.9));  // Bound within reasonable range

    // For win/loss ratio, use take profit vs stop loss
    double winAmount = positionSizing.takeProfitPct;
    double lossAmount = positionSizing.stopLossPct;

    if (signal.takeProfitPrice > 0 && signal.stopLossPrice > 0) {
        double currentPrice = ibConnector->getLatestPrice(signal.symbol).close;
        winAmount = std::abs(signal.takeProfitPrice - currentPrice) / currentPrice;
        lossAmount = std::abs(signal.stopLossPrice - currentPrice) / currentPrice;
    }

    // Ensure loss amount is positive
    lossAmount = std::max(lossAmount, 0.001);

    // Calculate win/loss ratio
    double winLossRatio = winAmount / lossAmount;

    // Kelly formula
    double kellyFraction = (winProbability * winLossRatio - (1.0 - winProbability)) / winLossRatio;

    // Bound kelly between 0 and max position size (conservative)
    return std::max(0.0, std::min(kellyFraction, positionSizing.maxPositionSize));
}

double TradeExecution::calculateValueAtRisk(
    const TradeSignal& signal, int timeHorizon, double confidenceLevel) const {

    // Calculate VaR using parametric method (assuming normal distribution)
    double volatility = signal.volatility;
    if (volatility <= 0.0) {
        // If volatility not provided, calculate from historical data
        PriceTimeSeries priceData = ibConnector->getHistoricalData(
            signal.symbol, "TODAY", 30, "1 day", "ADJUSTED_LAST"
        );

        std::vector<double> prices;
        for (const auto& data : priceData) {
            prices.push_back(data.close);
        }

        volatility = statistics->calculateVolatility(prices);
    }

    // Convert annual volatility to time horizon
    volatility = volatility / std::sqrt(252) * std::sqrt(timeHorizon);

    // Z-score for the confidence level (approximate for 95% = 1.645, 99% = 2.326)
    double zScore = 1.645;
    if (confidenceLevel > 0.95) {
        zScore = 2.326;
    } else if (confidenceLevel < 0.95) {
        zScore = 1.282;  // 90% confidence
    }

    // VaR as percentage of position
    double var = zScore * volatility;

    return var;
}

int TradeExecution::executeTrade(const TradeSignal& signal) {
    if (!ibConnector->isConnected()) {
        std::cerr << "Error: Not connected to Interactive Brokers" << std::endl;
        return -1;
    }

    // Check portfolio heat
    if (calculatePortfolioHeat() >= positionSizing.portfolioHeatMax) {
        std::cout << "Warning: Maximum portfolio heat reached, cannot execute trade" << std::endl;
        return -1;
    }

    // Check if we already have a position in this symbol
    auto it = currentPositions.find(signal.symbol);
    if (it != currentPositions.end() && it->second != 0) {
        std::cout << "Warning: Position already exists for " << signal.symbol << std::endl;
        return -1;
    }

    // Check if we've reached max positions
    if (currentPositions.size() >= positionSizing.maxPositions) {
        std::cout << "Warning: Maximum number of positions reached" << std::endl;
        return -1;
    }

    // Calculate position size
    int quantity = calculateShares(signal);
    if (quantity <= 0) {
        std::cout << "Warning: Calculated position size is zero or negative" << std::endl;
        return -1;
    }

    // Create order
    TradeOrder order;
    order.symbol = signal.symbol;
    order.quantity = quantity;
    order.price = ibConnector->getLatestPrice(signal.symbol).close;
    order.isBuy = signal.isBuy;
    order.orderType = "LIMIT";
    order.timeInForce = "DAY";

    // Place order
    int orderId = ibConnector->placeOrder(order);

    if (orderId > 0) {
        // Update position tracking (assume immediate execution for simplicity)
        currentPositions[signal.symbol] = signal.isBuy ? quantity : -quantity;
        positionValues[signal.symbol] = quantity * order.price;
        availableCapital -= quantity * order.price;

        // Update portfolio metrics
        updatePortfolioMetrics();

        std::cout << "Trade executed: " << (signal.isBuy ? "BUY " : "SELL ")
                  << quantity << " shares of " << signal.symbol
                  << " at $" << order.price << std::endl;
    }

    return orderId;
}

std::pair<int, int> TradeExecution::executePairTrade(
    const TradeSignal& signal, const TradeSignal& hedgeSignal) {

    if (!ibConnector->isConnected()) {
        std::cerr << "Error: Not connected to Interactive Brokers" << std::endl;
        return {-1, -1};
    }

    // Check portfolio heat
    if (calculatePortfolioHeat() >= positionSizing.portfolioHeatMax) {
        std::cout << "Warning: Maximum portfolio heat reached, cannot execute pair trade" << std::endl;
        return {-1, -1};
    }

    // Calculate position sizes
    auto [sharesMain, sharesHedge] = calculatePairShares(signal, hedgeSignal);

    if (sharesMain <= 0 || sharesHedge <= 0) {
        std::cout << "Warning: Calculated position sizes are zero or negative" << std::endl;
        return {-1, -1};
    }

    // Create main leg order
    TradeOrder mainOrder;
    mainOrder.symbol = signal.symbol;
    mainOrder.quantity = sharesMain;
    mainOrder.price = ibConnector->getLatestPrice(signal.symbol).close;
    mainOrder.isBuy = signal.isBuy;
    mainOrder.orderType = "LIMIT";
    mainOrder.timeInForce = "DAY";

    // Create hedge leg order
    TradeOrder hedgeOrder;
    hedgeOrder.symbol = hedgeSignal.symbol;
    hedgeOrder.quantity = sharesHedge;
    hedgeOrder.price = ibConnector->getLatestPrice(hedgeSignal.symbol).close;
    hedgeOrder.isBuy = !signal.isBuy;  // Opposite direction for hedge
    hedgeOrder.orderType = "LIMIT";
    hedgeOrder.timeInForce = "DAY";

    // Place orders
    int mainOrderId = ibConnector->placeOrder(mainOrder);
    int hedgeOrderId = ibConnector->placeOrder(hedgeOrder);

    if (mainOrderId > 0 && hedgeOrderId > 0) {
        // Update position tracking (assume immediate execution for simplicity)
        int mainDirection = signal.isBuy ? 1 : -1;
        int hedgeDirection = !signal.isBuy ? 1 : -1;

        currentPositions[signal.symbol] = sharesMain * mainDirection;
        currentPositions[hedgeSignal.symbol] = sharesHedge * hedgeDirection;

        positionValues[signal.symbol] = sharesMain * mainOrder.price;
        positionValues[hedgeSignal.symbol] = sharesHedge * hedgeOrder.price;

        availableCapital -= (sharesMain * mainOrder.price + sharesHedge * hedgeOrder.price);

        // Update portfolio metrics
        updatePortfolioMetrics();

        std::cout << "Pair trade executed:" << std::endl;
        std::cout << "  " << (signal.isBuy ? "LONG " : "SHORT ")
                  << sharesMain << " shares of " << signal.symbol
                  << " at $" << mainOrder.price << std::endl;
        std::cout << "  " << (!signal.isBuy ? "LONG " : "SHORT ")
                  << sharesHedge << " shares of " << hedgeSignal.symbol
                  << " at $" << hedgeOrder.price << std::endl;
    }

    return {mainOrderId, hedgeOrderId};
}

bool TradeExecution::closePosition(const std::string& symbol) {
    auto it = currentPositions.find(symbol);
    if (it == currentPositions.end() || it->second == 0) {
        std::cout << "No position found for " << symbol << std::endl;
        return false;
    }

    // Create order to close position
    TradeOrder order;
    order.symbol = symbol;
    order.quantity = std::abs(it->second);
    order.price = ibConnector->getLatestPrice(symbol).close;
    order.isBuy = it->second < 0;  // Buy if short, sell if long
    order.orderType = "MARKET";
    order.timeInForce = "DAY";

    // Place order
    int orderId = ibConnector->placeOrder(order);

    if (orderId > 0) {
        // Update tracking (assume immediate execution)
        availableCapital += positionValues[symbol];
        currentPositions[symbol] = 0;
        positionValues[symbol] = 0.0;

        // Update portfolio metrics
        updatePortfolioMetrics();

        std::cout << "Position closed: " << symbol << std::endl;
        return true;
    }

    return false;
}

bool TradeExecution::closeAllPositions() {
    bool success = true;

    std::vector<std::string> openPositions;
    for (const auto& [symbol, quantity] : currentPositions) {
        if (quantity != 0) {
            openPositions.push_back(symbol);
        }
    }

    for (const auto& symbol : openPositions) {
        if (!closePosition(symbol)) {
            success = false;
        }
    }

    return success;
}

void TradeExecution::setPositionSizingModel(PositionSizingModel model) {
    positionSizing.model = model;
    std::cout << "Position sizing model updated" << std::endl;
}

void TradeExecution::setMaxPositionSize(double maxPositionSize) {
    positionSizing.maxPositionSize = maxPositionSize;
    std::cout << "Max position size updated to " << maxPositionSize * 100.0 << "%" << std::endl;
}

void TradeExecution::setRiskPerTrade(double riskPerTrade) {
    positionSizing.riskPerTrade = riskPerTrade;
    std::cout << "Risk per trade updated to " << riskPerTrade * 100.0 << "%" << std::endl;
}
