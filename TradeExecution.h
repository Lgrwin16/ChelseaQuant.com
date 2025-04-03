#ifndef TRADE_EXECUTION_H
#define TRADE_EXECUTION_H

#include <string>
#include <vector>
#include <map>
#include "IBConnector.h"
#include "Statistics.h"

enum class PositionSizingModel {
    FIXED_PERCENTAGE,     // Fixed percentage of capital
    VOLATILITY_ADJUSTED,  // Adjusts size based on volatility
    KELLY_CRITERION,      // Optimal f (Kelly criterion)
    EQUAL_RISK,           // Equal risk across all positions
    VALUE_AT_RISK         // VaR-based position sizing
};

struct PositionSizing {
    PositionSizingModel model;
    double maxPositionSize;        // Maximum position size as percentage of capital (e.g., 0.05 for 5%)
    double riskPerTrade;           // Maximum risk per trade as percentage of capital
    double portfolioHeatMax;       // Maximum portfolio heat (total risk exposure) as percentage
    int maxPositions;              // Maximum number of open positions
    double stopLossPct;            // Stop loss percentage
    double takeProfitPct;          // Take profit percentage

    // Default constructor with conservative settings
    PositionSizing()
        : model(PositionSizingModel::VOLATILITY_ADJUSTED),
          maxPositionSize(0.05),
          riskPerTrade(0.01),
          portfolioHeatMax(0.20),
          maxPositions(10),
          stopLossPct(0.02),
          takeProfitPct(0.03) {}
};

struct TradeSignal {
    std::string symbol;
    bool isBuy;
    double targetPrice;
    double confidenceScore;
    double volatility;
    double beta;
    double stopLossPrice;
    double takeProfitPrice;

    // Optional fields for pair trading
    std::string pairedSymbol;
    double hedgeRatio;
};

class TradeExecution {
private:
    std::shared_ptr<IBConnector> ibConnector;
    std::shared_ptr<Statistics> statistics;
    PositionSizing positionSizing;
    double accountCapital;
    double availableCapital;

    // Portfolio management
    std::map<std::string, int> currentPositions;
    std::map<std::string, double> positionValues;
    double totalPortfolioValue;
    double totalRiskExposure;

    // Private methods
    double calculatePositionSize(const TradeSignal& signal) const;
    double calculateKellyCriterion(const TradeSignal& signal) const;
    double calculateVolatilityAdjustedSize(const TradeSignal& signal) const;
    double calculateValueAtRisk(const TradeSignal& signal, int timeHorizon, double confidenceLevel) const;
    void updatePortfolioMetrics();

public:
    TradeExecution(
        std::shared_ptr<IBConnector> ibConnector,
        std::shared_ptr<Statistics> statistics,
        double initialCapital = 1000000.0,
        const PositionSizing& positionSizing = PositionSizing()
    );

    ~TradeExecution();

    // Portfolio management
    void updateAccountInfo();
    double getAccountCapital() const;
    double getAvailableCapital() const;
    double getTotalPortfolioValue() const;
    double getTotalRiskExposure() const;

    // Position sizing
    int calculateShares(const TradeSignal& signal);
    std::pair<int, int> calculatePairShares(const TradeSignal& signal, const TradeSignal& hedgeSignal);

    // Order execution
    int executeTrade(const TradeSignal& signal);
    std::pair<int, int> executePairTrade(const TradeSignal& signal, const TradeSignal& hedgeSignal);
    bool closePosition(const std::string& symbol);
    bool closeAllPositions();

    // Risk management
    void setPositionSizingModel(PositionSizingModel model);
    void setMaxPositionSize(double maxPositionSize);
    void setRiskPerTrade(double riskPerTrade);
    double calculatePortfolioHeat() const;
};

#endif // TRADE_EXECUTION_H