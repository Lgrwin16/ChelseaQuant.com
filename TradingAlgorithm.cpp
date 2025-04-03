#include "TradingAlgorithm.h"
#include <iostream>
#include <algorithm>
#include <cmath>

TradingAlgorithm::TradingAlgorithm(TradingStrategy strategy,
                                 const std::vector<std::string>& symbols,
                                 const std::string& marketIndex,
                                 double riskFreeRate,
                                 double lookbackPeriod,
                                 double confidenceThreshold)
    : strategy(strategy),
      symbols(symbols),
      marketIndex(marketIndex),
      riskFreeRate(riskFreeRate),
      lookbackPeriod(lookbackPeriod),
      confidenceThreshold(confidenceThreshold) {

}

TradingAlgorithm::~TradingAlgorithm() {
    // Clean up resources if needed
}

void TradingAlgorithm::initialize() {
    // Initialize IB connector
    ibConnector = std::make_shared<IBConnector>();
    if (!ibConnector->connect()) {
        std::cerr << "Failed to connect to Interactive Brokers" << std::endl;
        return;
    }

    // Initialize market data
    marketData = std::make_shared<MarketData>(marketIndex, static_cast<int>(lookbackPeriod));

    // Initialize statistics engine
    statistics = std::make_shared<Statistics>(marketData);

    std::cout << "Trading algorithm initialized" << std::endl;
}

void TradingAlgorithm::fetchData() {
    // Add market index to the list of symbols to fetch
    std::vector<std::string> allSymbols = symbols;
    if (std::find(allSymbols.begin(), allSymbols.end(), marketIndex) == allSymbols.end()) {
        allSymbols.push_back(marketIndex);
    }

    // Fetch historical data
    marketData->fetchHistoricalData(allSymbols);

    // Update with the latest live data
    marketData->updateLiveData(allSymbols);

    std::cout << "Market data fetched for " << allSymbols.size() << " symbols" << std::endl;
}

void TradingAlgorithm::analyzeData() {
    // Calculate CAPM parameters (betas and alphas)
    calculateBetasAndAlphas();

    // Calculate mean-reversion parameters
    calculateMeanReversionParameters();

    // Calculate stochastic process parameters
    calculateStochasticParameters();

    std::cout << "Data analysis completed" << std::endl;
}

void TradingAlgorithm::calculateBetasAndAlphas() {
    for (const auto& symbol : symbols) {
        // Calculate beta using CAPM
        double beta = statistics->calculateBeta(symbol);
        betas[symbol] = beta;

        // Calculate alpha (excess return)
        double alpha = statistics->calculateAlpha(symbol, riskFreeRate);
        alphas[symbol] = alpha;

        std::cout << "Symbol: " << symbol << ", Beta: " << beta << ", Alpha: " << alpha << std::endl;
    }
}

void TradingAlgorithm::calculateMeanReversionParameters() {
    for (const auto& symbol : symbols) {
        PriceTimeSeries priceData = marketData->getStockData(symbol);
        std::vector<double> prices;

        for (const auto& data : priceData) {
            prices.push_back(data.close);
        }

        // Calculate mean price
        meanPrices[symbol] = statistics->mean(prices);

        // Calculate standard deviation
        stdDeviations[symbol] = statistics->standardDeviation(prices);

        std::cout << "Symbol: " << symbol
                  << ", Mean: " << meanPrices[symbol]
                  << ", StdDev: " << stdDeviations[symbol] << std::endl;
    }
}

void TradingAlgorithm::calculateStochasticParameters() {
    for (const auto& symbol : symbols) {
        // Get price data
        std::vector<double> returns = marketData->getReturns(symbol);

        // Calculate drift (μ) and volatility (σ) for GBM
        drifts[symbol] = statistics->calculateDrift(returns);
        volatilities[symbol] = statistics->calculateVolatility(returns);

        std::cout << "Symbol: " << symbol
                  << ", Drift: " << drifts[symbol]
                  << ", Volatility: " << volatilities[symbol] << std::endl;
    }
}

void TradingAlgorithm::generateTrades() {
    if (strategy == TradingStrategy::STAT_ARB) {
        runStatArbStrategy();
    } else {
        runMeanReversionStrategy();
    }
}

void TradingAlgorithm::runStatArbStrategy() {
    std::cout << "Running Statistical Arbitrage Strategy..." << std::endl;

    // Sort stocks by alpha (descending)
    std::vector<std::pair<std::string, double>> sortedAlphas;
    for (const auto& alphaPair : alphas) {
        sortedAlphas.push_back(alphaPair);
    }

    std::sort(sortedAlphas.begin(), sortedAlphas.end(),
             [](const auto& a, const auto& b) { return a.second > b.second; });

    // Select top stocks with positive alpha and bottom stocks with negative alpha
    std::vector<std::string> longCandidates;
    std::vector<std::string> shortCandidates;

    for (const auto& [symbol, alpha] : sortedAlphas) {
        // Calculate Sharpe ratio for additional filtering
        double sharpe = calculateSharpeRatio(symbol);
        double infoRatio = calculateInformationRatio(symbol);

        if (alpha > confidenceThreshold && sharpe > 1.0 && infoRatio > 0.5) {
            longCandidates.push_back(symbol);
            std::cout << "Long candidate: " << symbol
                      << ", Alpha: " << alpha
                      << ", Sharpe: " << sharpe
                      << ", Info Ratio: " << infoRatio << std::endl;
        } else if (alpha < -confidenceThreshold && sharpe < -0.5 && infoRatio < -0.3) {
            shortCandidates.push_back(symbol);
            std::cout << "Short candidate: " << symbol
                      << ", Alpha: " << alpha
                      << ", Sharpe: " << sharpe
                      << ", Info Ratio: " << infoRatio << std::endl;
        }
    }

    // Generate trade orders
    for (const auto& symbol : longCandidates) {
        // Calculate position size based on volatility or other metrics
        double currentPrice = marketData->getCurrentPrice(symbol);
        int quantity = static_cast<int>(10000 / currentPrice); // Simple position sizing

        TradeOrder order;
        order.symbol = symbol;
        order.quantity = quantity;
        order.price = currentPrice;
        order.isBuy = true;
        order.orderType = "LIMIT";
        order.timeInForce = "DAY";

        int orderId = ibConnector->placeOrder(order);
        std::cout << "Placed BUY order for " << symbol
                  << ", Quantity: " << quantity
                  << ", Price: " << currentPrice
                  << ", OrderID: " << orderId << std::endl;
    }

    for (const auto& symbol : shortCandidates) {
        // Calculate position size
        double currentPrice = marketData->getCurrentPrice(symbol);
        int quantity = static_cast<int>(5000 / currentPrice); // Less exposure for shorts

        TradeOrder order;
        order.symbol = symbol;
        order.quantity = quantity;
        order.price = currentPrice;
        order.isBuy = false;
        order.orderType = "LIMIT";
        order.timeInForce = "DAY";

        int orderId = ibConnector->placeOrder(order);
        std::cout << "Placed SELL order for " << symbol
                  << ", Quantity: " << quantity
                  << ", Price: " << currentPrice
                  << ", OrderID: " << orderId << std::endl;
    }
}

void TradingAlgorithm::runMeanReversionStrategy() {
    std::cout << "Running Mean-Reversion Strategy..." << std::endl;

    for (const auto& symbol : symbols) {
        double currentPrice = marketData->getCurrentPrice(symbol);
        double mean = meanPrices[symbol];
        double stdDev = stdDeviations[symbol];

        // Calculate z-score
        double zScore = (currentPrice - mean) / stdDev;

        // Run stationarity test
        std::vector<double> prices;
        PriceTimeSeries priceData = marketData->getStockData(symbol);
        for (const auto& data : priceData) {
            prices.push_back(data.close);
        }

        double adfStat = statistics->augmentedDickeyFullerTest(prices);
        double hurstExponent = statistics->hurst(prices);
        double halfLifeValue = statistics->halfLife(prices);

        std::cout << "Symbol: " << symbol
                  << ", Z-Score: " << zScore
                  << ", ADF: " << adfStat
                  << ", Hurst: " << hurstExponent
                  << ", Half-life: " << halfLifeValue << std::endl;

        // Trading logic for mean-reversion
        bool isStationary = adfStat < -2.9; // ADF critical value
        bool isMeanReverting = hurstExponent < 0.5 && halfLifeValue > 0 && halfLifeValue < 20;

        if (isStationary && isMeanReverting) {
            TradeOrder order;
            order.symbol = symbol;
            int quantity = static_cast<int>(10000 / currentPrice);
            order.quantity = quantity;
            order.price = currentPrice;
            order.timeInForce = "DAY";

            if (zScore > confidenceThreshold) {
                // Significantly above mean - go short
                order.isBuy = false;
                order.orderType = "LIMIT";
                int orderId = ibConnector->placeOrder(order);
                std::cout << "Placed SELL order for " << symbol
                          << ", Quantity: " << quantity
                          << ", Price: " << currentPrice
                          << ", OrderID: " << orderId << std::endl;
            } else if (zScore < -confidenceThreshold) {
                // Significantly below mean - go long
                order.isBuy = true;
                order.orderType = "LIMIT";
                int orderId = ibConnector->placeOrder(order);
                std::cout << "Placed BUY order for " << symbol
                          << ", Quantity: " << quantity
                          << ", Price: " << currentPrice
                          << ", OrderID: " << orderId << std::endl;
            }
        }
    }
}

void TradingAlgorithm::executeTrades() {
    // This would monitor open orders and handle executions
    // For simplicity, we'll just check the status of pending orders
    std::map<int, TradeOrder> openOrders = ibConnector->getOpenOrders();

    std::cout << "Monitoring " << openOrders.size() << " open orders" << std::endl;

    // In a real system, we would implement callbacks for order status updates
}

double TradingAlgorithm::calculateSharpeRatio(const std::string& symbol) {
    std::vector<double> returns = marketData->getReturns(symbol);
    return statistics->sharpeRatio(returns, riskFreeRate);
}

double TradingAlgorithm::calculateInformationRatio(const std::string& symbol) {
    std::vector<double> returns = marketData->getReturns(symbol);
    std::vector<double> benchmarkReturns = marketData->getMarketReturns();
    return statistics->informationRatio(returns, benchmarkReturns);
}

double TradingAlgorithm::expectedReturn(const std::string& symbol) {
    // Using CAPM to calculate expected return
    std::vector<double> marketReturns = marketData->getMarketReturns();
    double avgMarketReturn = statistics->mean(marketReturns);

    return statistics->expectedReturn(symbol, avgMarketReturn, riskFreeRate);
}

void TradingAlgorithm::runBacktest(const std::string& startDate, const std::string& endDate) {
    std::cout << "Running backtest from " << startDate << " to " << endDate << std::endl;
    // Backtest implementation would go here
    // This would involve historical data simulation and performance analysis
}

std::map<std::string, double> TradingAlgorithm::getAlphas() const {
    return alphas;
}

std::map<std::string, double> TradingAlgorithm::getBetas() const {
    return betas;
}

TradingStrategy TradingAlgorithm::getStrategy() const {
    return strategy;
}