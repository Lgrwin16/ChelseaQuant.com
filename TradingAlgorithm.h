#ifndef TRADING_ALGORITHM_H
#define TRADING_ALGORITHM_H

#include <vector>
#include <string>
#include <map>
#include <memory>
#include "MarketData.h"
#include "Statistics.h"
#include "IBConnector.h"

enum class TradingStrategy {
    STAT_ARB,
    MEAN_REVERSION
};

class TradingAlgorithm {
private:
    // Configuration
    TradingStrategy strategy;
    std::vector<std::string> symbols;
    double riskFreeRate;
    double lookbackPeriod;
    double confidenceThreshold;

    // Market data and results
    std::shared_ptr<IBConnector> ibConnector;
    std::shared_ptr<MarketData> marketData;
    std::shared_ptr<Statistics> statistics;

    // CAPM model parameters
    std::map<std::string, double> betas;
    std::map<std::string, double> alphas;
    std::string marketIndex;

    // Mean-reversion parameters
    std::map<std::string, double> meanPrices;
    std::map<std::string, double> stdDeviations;

    // Stochastic process parameters
    std::map<std::string, double> drifts;
    std::map<std::string, double> volatilities;

    // Private methods
    void calculateBetasAndAlphas();
    void calculateMeanReversionParameters();
    void calculateStochasticParameters();
    void runStatArbStrategy();
    void runMeanReversionStrategy();
    double calculateSharpeRatio(const std::string& symbol);
    double calculateInformationRatio(const std::string& symbol);
    double expectedReturn(const std::string& symbol);

public:
    TradingAlgorithm(TradingStrategy strategy,
                    const std::vector<std::string>& symbols,
                    const std::string& marketIndex,
                    double riskFreeRate,
                    double lookbackPeriod,
                    double confidenceThreshold);

    ~TradingAlgorithm();

    void initialize();
    void fetchData();
    void analyzeData();
    void generateTrades();
    void executeTrades();
    void runBacktest(const std::string& startDate, const std::string& endDate);

    // Getters
    std::map<std::string, double> getAlphas() const;
    std::map<std::string, double> getBetas() const;
    TradingStrategy getStrategy() const;
};

#endif // TRADING_ALGORITHM_H