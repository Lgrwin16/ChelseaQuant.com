#ifndef STATISTICS_H
#define STATISTICS_H

#include <vector>
#include <string>
#include <map>
#include <random>
#include "MarketData.h"

class Statistics {
private:
    std::shared_ptr<MarketData> marketData;
    std::mt19937_64 rng; // Random number generator for stochastic processes

public:
    Statistics(std::shared_ptr<MarketData> marketData);
    ~Statistics();

    // Basic statistics
    double mean(const std::vector<double>& data) const;
    double standardDeviation(const std::vector<double>& data) const;
    double covariance(const std::vector<double>& x, const std::vector<double>& y) const;
    double correlation(const std::vector<double>& x, const std::vector<double>& y) const;

    // CAPM model
    double calculateBeta(const std::string& symbol) const;
    double calculateAlpha(const std::string& symbol, double riskFreeRate) const;
    double expectedReturn(const std::string& symbol, double marketReturn, double riskFreeRate) const;

    // Mean-reversion tests
    double halfLife(const std::vector<double>& data) const;
    double augmentedDickeyFullerTest(const std::vector<double>& data) const;
    double hurst(const std::vector<double>& data) const;

    // Stochastic calculus
    std::vector<double> geometricBrownianMotion(double s0, double mu, double sigma, int steps, double dt) const;
    double calculateDrift(const std::vector<double>& prices) const;
    double calculateVolatility(const std::vector<double>& prices) const;
    double calculateOptionPrice(double s, double k, double t, double r, double sigma, bool isCall) const;

    // Portfolio metrics
    double sharpeRatio(const std::vector<double>& returns, double riskFreeRate) const;
    double informationRatio(const std::vector<double>& returns, const std::vector<double>& benchmark) const;
    double maxDrawdown(const std::vector<double>& equityCurve) const;

    // Pair trading metrics
    double zScore(const std::vector<double>& spread) const;
    std::vector<double> calculateSpread(const std::string& symbol1, const std::string& symbol2, double hedge) const;
    double optimalHedgeRatio(const std::string& symbol1, const std::string& symbol2) const;

    // Monte Carlo simulations
    std::vector<std::vector<double>> monteCarloSimulation(const std::string& symbol, int numSimulations, int timeSteps) const;
};

#endif // STATISTICS_H