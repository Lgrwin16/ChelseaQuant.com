#include "Statistics.h"
#include <cmath>
#include <numeric>
#include <algorithm>
#include <random>
#include <iostream>

Statistics::Statistics(std::shared_ptr<MarketData> marketData)
    : marketData(marketData), rng(std::random_device{}()) {
}

Statistics::~Statistics() {
    // Clean up resources if needed
}

double Statistics::mean(const std::vector<double>& data) const {
    if (data.empty()) {
        return 0.0;
    }
    return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
}

double Statistics::standardDeviation(const std::vector<double>& data) const {
    if (data.size() <= 1) {
        return 0.0;
    }

    double avg = mean(data);
    double sum = 0.0;

    for (const auto& value : data) {
        sum += (value - avg) * (value - avg);
    }

    return std::sqrt(sum / (data.size() - 1));
}

double Statistics::covariance(const std::vector<double>& x, const std::vector<double>& y) const {
    if (x.size() != y.size() || x.empty()) {
        return 0.0;
    }

    double meanX = mean(x);
    double meanY = mean(y);
    double sum = 0.0;

    for (size_t i = 0; i < x.size(); ++i) {
        sum += (x[i] - meanX) * (y[i] - meanY);
    }

    return sum / (x.size() - 1);
}

double Statistics::correlation(const std::vector<double>& x, const std::vector<double>& y) const {
    double sdX = standardDeviation(x);
    double sdY = standardDeviation(y);

    if (sdX == 0.0 || sdY == 0.0) {
        return 0.0;
    }

    return covariance(x, y) / (sdX * sdY);
}

double Statistics::calculateBeta(const std::string& symbol) const {
    std::vector<double> stockReturns = marketData->getReturns(symbol);
    std::vector<double> marketReturns = marketData->getMarketReturns();

    // Ensure same length for both return series
    size_t minSize = std::min(stockReturns.size(), marketReturns.size());
    stockReturns.resize(minSize);
    marketReturns.resize(minSize);

    // Calculate beta = Cov(stock, market) / Var(market)
    double covar = covariance(stockReturns, marketReturns);
    double marketVar = covariance(marketReturns, marketReturns);

    if (marketVar == 0.0) {
        return 0.0;
    }

    return covar / marketVar;
}

double Statistics::calculateAlpha(const std::string& symbol, double riskFreeRate) const {
    std::vector<double> stockReturns = marketData->getReturns(symbol);
    std::vector<double> marketReturns = marketData->getMarketReturns();

    if (stockReturns.empty() || marketReturns.empty()) {
        return 0.0;
    }

    // Calculate average returns
    double avgStockReturn = mean(stockReturns);
    double avgMarketReturn = mean(marketReturns);

    // Calculate beta
    double beta = calculateBeta(symbol);

    // Calculate annualized returns (assuming daily returns)
    double annualizedStockReturn = avgStockReturn * 252;
    double annualizedMarketReturn = avgMarketReturn * 252;

    // Calculate alpha using CAPM:
    // Alpha = R_stock - [R_f + Beta * (R_market - R_f)]
    double alpha = annualizedStockReturn - (riskFreeRate + beta * (annualizedMarketReturn - riskFreeRate));

    return alpha;
}

double Statistics::expectedReturn(const std::string& symbol, double marketReturn, double riskFreeRate) const {
    // CAPM expected return formula: E(R) = R_f + β(R_m - R_f)
    double beta = calculateBeta(symbol);
    return riskFreeRate + beta * (marketReturn - riskFreeRate);
}

double Statistics::halfLife(const std::vector<double>& data) const {
    if (data.size() < 3) {
        return 0.0;
    }

    // Create vectors for regression: y = price changes, x = lagged prices
    std::vector<double> y, x;
    for (size_t i = 1; i < data.size(); ++i) {
        y.push_back(data[i] - data[i-1]);  // Price change
        x.push_back(data[i-1]);            // Lagged price
    }

    // Calculate regression coefficient (phi)
    double sumX = std::accumulate(x.begin(), x.end(), 0.0);
    double sumY = std::accumulate(y.begin(), y.end(), 0.0);
    double sumXY = 0.0, sumX2 = 0.0;

    for (size_t i = 0; i < x.size(); ++i) {
        sumXY += x[i] * y[i];
        sumX2 += x[i] * x[i];
    }

    double meanX = sumX / x.size();
    double meanY = sumY / y.size();
    double phi = (sumXY - x.size() * meanX * meanY) / (sumX2 - x.size() * meanX * meanX);

    // Half-life = ln(0.5) / ln(1 + phi)
    if (phi >= 0) {
        return 0.0;  // Not mean-reverting
    }

    return std::log(0.5) / std::log(1 + phi);
}

double Statistics::augmentedDickeyFullerTest(const std::vector<double>& data) const {
    // This is a simplified ADF test implementation
    // In a production system, you'd use a statistical library

    if (data.size() < 20) {  // Need sufficient data for meaningful test
        return 0.0;
    }

    // Create vectors for regression: y = price changes, x = lagged prices
    std::vector<double> y, x;
    for (size_t i = 1; i < data.size(); ++i) {
        y.push_back(data[i] - data[i-1]);  // Price change
        x.push_back(data[i-1]);            // Lagged price
    }

    // Calculate regression coefficient (phi)
    double sumX = std::accumulate(x.begin(), x.end(), 0.0);
    double sumY = std::accumulate(y.begin(), y.end(), 0.0);
    double sumXY = 0.0, sumX2 = 0.0;

    for (size_t i = 0; i < x.size(); ++i) {
        sumXY += x[i] * y[i];
        sumX2 += x[i] * x[i];
    }

    double meanX = sumX / x.size();
    double meanY = sumY / y.size();
    double phi = (sumXY - x.size() * meanX * meanY) / (sumX2 - x.size() * meanX * meanX);

    // Calculate standard error (simplified)
    double sse = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        double predicted = meanY + phi * (x[i] - meanX);
        double error = y[i] - predicted;
        sse += error * error;
    }

    double variance = sse / (x.size() - 2);
    double stdError = std::sqrt(variance / sumX2);

    // Calculate ADF statistic
    double adfStat = phi / stdError;

    return adfStat;
}

double Statistics::hurst(const std::vector<double>& data) const {
    if (data.size() < 20) {
        return 0.5;  // Default for insufficient data
    }

    const int maxLag = std::min(static_cast<int>(data.size() / 4), 20);
    std::vector<double> logLags;
    std::vector<double> logRSs;

    for (int lag = 10; lag <= maxLag; ++lag) {
        std::vector<double> rs;

        // Calculate R/S for different segments
        for (size_t i = 0; i <= data.size() - lag; i += lag / 2) {
            if (i + lag > data.size()) continue;

            // Extract segment
            std::vector<double> segment(data.begin() + i, data.begin() + i + lag);

            // Calculate mean
            double segmentMean = mean(segment);

            // Calculate cumulative deviations
            std::vector<double> cumDevs(lag, 0.0);
            for (int j = 0; j < lag; ++j) {
                cumDevs[j] = segment[j] - segmentMean;
                if (j > 0) cumDevs[j] += cumDevs[j-1];
            }

            // Find range (R)
            double maxDev = *std::max_element(cumDevs.begin(), cumDevs.end());
            double minDev = *std::min_element(cumDevs.begin(), cumDevs.end());
            double range = maxDev - minDev;

            // Calculate standard deviation (S)
            double stdDev = standardDeviation(segment);

            // Calculate R/S ratio
            if (stdDev > 0) {
                rs.push_back(range / stdDev);
            }
        }

        // Calculate average R/S for this lag
        if (!rs.empty()) {
            double avgRS = mean(rs);
            logLags.push_back(std::log(lag));
            logRSs.push_back(std::log(avgRS));
        }
    }

    // Linear regression on log-log plot to estimate Hurst exponent
    if (logLags.size() < 4) {
        return 0.5;  // Default if not enough points for regression
    }

    // Calculate slope (Hurst exponent)
    double sumX = std::accumulate(logLags.begin(), logLags.end(), 0.0);
    double sumY = std::accumulate(logRSs.begin(), logRSs.end(), 0.0);
    double sumXY = 0.0, sumX2 = 0.0;

    for (size_t i = 0; i < logLags.size(); ++i) {
        sumXY += logLags[i] * logRSs[i];
        sumX2 += logLags[i] * logLags[i];
    }

    int n = logLags.size();
    double hurst = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);

    return hurst;
}

std::vector<double> Statistics::geometricBrownianMotion(
    double s0, double mu, double sigma, int steps, double dt) const {

    std::vector<double> path(steps + 1);
    path[0] = s0;

    std::normal_distribution<> normalDist(0.0, 1.0);

    for (int i = 0; i < steps; ++i) {
        double z = normalDist(rng);
        path[i+1] = path[i] * std::exp((mu - 0.5 * sigma * sigma) * dt + sigma * std::sqrt(dt) * z);
    }

    return path;
}

double Statistics::calculateDrift(const std::vector<double>& prices) const {
    if (prices.size() < 2) {
        return 0.0;
    }

    // Calculate log returns
    std::vector<double> logReturns;
    for (size_t i = 1; i < prices.size(); ++i) {
        logReturns.push_back(std::log(prices[i-1] / prices[i]));
    }

    // Drift is the mean of log returns (adjusted for time and volatility)
    double meanReturn = mean(logReturns);
    double vol = standardDeviation(logReturns);

    // Annualize (assuming daily data)
    return meanReturn * 252 + 0.5 * vol * vol;
}

double Statistics::calculateVolatility(const std::vector<double>& prices) const {
    if (prices.size() < 2) {
        return 0.0;
    }

    // Calculate log returns
    std::vector<double> logReturns;
    for (size_t i = 1; i < prices.size(); ++i) {
        logReturns.push_back(std::log(prices[i-1] / prices[i]));
    }

    // Volatility is the standard deviation of log returns
    double vol = standardDeviation(logReturns);

    // Annualize (assuming daily data)
    return vol * std::sqrt(252);
}

double Statistics::calculateOptionPrice(
    double s, double k, double t, double r, double sigma, bool isCall) const {

    // Black-Scholes option pricing model
    double d1 = (std::log(s / k) + (r + sigma * sigma / 2) * t) / (sigma * std::sqrt(t));
    double d2 = d1 - sigma * std::sqrt(t);

    // Normal CDF approximation (Abramowitz and Stegun method)
    auto normalCDF = [](double x) {
        const double a1 = 0.254829592;
        const double a2 = -0.284496736;
        const double a3 = 1.421413741;
        const double a4 = -1.453152027;
        const double a5 = 1.061405429;
        const double p = 0.3275911;

        int sign = 1;
        if (x < 0) {
            sign = -1;
            x = -x;
        }

        double t = 1.0 / (1.0 + p * x);
        double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(-x * x);

        return 0.5 * (1.0 + sign * y);
    };

    if (isCall) {
        return s * normalCDF(d1) - k * std::exp(-r * t) * normalCDF(d2);
    } else {
        return k * std::exp(-r * t) * normalCDF(-d2) - s * normalCDF(-d1);
    }
}

double Statistics::sharpeRatio(const std::vector<double>& returns, double riskFreeRate) const {
    if (returns.empty()) {
        return 0.0;
    }

    double avgReturn = mean(returns);
    double stdDev = standardDeviation(returns);

    if (stdDev == 0.0) {
        return 0.0;
    }

    // Calculate annualized Sharpe ratio (assuming daily returns)
    double excessReturn = avgReturn - riskFreeRate / 252;  // Daily risk-free rate
    double annualizedExcessReturn = excessReturn * 252;
    double annualizedStdDev = stdDev * std::sqrt(252);

    return annualizedExcessReturn / annualizedStdDev;
}

double Statistics::informationRatio(
    const std::vector<double>& returns, const std::vector<double>& benchmark) const {

    if (returns.empty() || benchmark.empty() || returns.size() != benchmark.size()) {
        return 0.0;
    }

    // Calculate excess returns
    std::vector<double> excessReturns;
    for (size_t i = 0; i < returns.size(); ++i) {
        excessReturns.push_back(returns[i] - benchmark[i]);
    }

    double avgExcessReturn = mean(excessReturns);
    double trackingError = standardDeviation(excessReturns);

    if (trackingError == 0.0) {
        return 0.0;
    }

    // Annualize (assuming daily returns)
    double annualizedExcessReturn = avgExcessReturn * 252;
    double annualizedTrackingError = trackingError * std::sqrt(252);

    return annualizedExcessReturn / annualizedTrackingError;
}

double Statistics::maxDrawdown(const std::vector<double>& equityCurve) const {
    if (equityCurve.size() < 2) {
        return 0.0;
    }

    double maxDrawdown = 0.0;
    double peak = equityCurve[0];

    for (size_t i = 1; i < equityCurve.size(); ++i) {
        if (equityCurve[i] > peak) {
            peak = equityCurve[i];
        } else {
            double drawdown = (peak - equityCurve[i]) / peak;
            if (drawdown > maxDrawdown) {
                maxDrawdown = drawdown;
            }
        }
    }

    return maxDrawdown;
}

double Statistics::zScore(const std::vector<double>& spread) const {
    if (spread.size() < 2) {
        return 0.0;
    }

    double spreadMean = mean(spread);
    double spreadStd = standardDeviation(spread);

    if (spreadStd == 0.0) {
        return 0.0;
    }

    return (spread.back() - spreadMean) / spreadStd;
}

std::vector<double> Statistics::calculateSpread(
    const std::string& symbol1, const std::string& symbol2, double hedge) const {

    PriceTimeSeries prices1 = marketData->getStockData(symbol1);
    PriceTimeSeries prices2 = marketData->getStockData(symbol2);

    // Ensure same length
    size_t minSize = std::min(prices1.size(), prices2.size());

    std::vector<double> spread;
    for (size_t i = 0; i < minSize; ++i) {
        spread.push_back(prices1[i].close - hedge * prices2[i].close);
    }

    return spread;
}

double Statistics::optimalHedgeRatio(const std::string& symbol1, const std::string& symbol2) const {
    std::vector<double> returns1 = marketData->getReturns(symbol1);
    std::vector<double> returns2 = marketData->getReturns(symbol2);

    // OLS regression to find optimal hedge ratio
    double beta = covariance(returns1, returns2) / covariance(returns2, returns2);

    return beta;
}

std::vector<std::vector<double>> Statistics::monteCarloSimulation(
    const std::string& symbol, int numSimulations, int timeSteps) const {

    double currentPrice = marketData->getCurrentPrice(symbol);
    double mu = calculateDrift(marketData->getStockData(symbol));
    double sigma = calculateVolatility(marketData->getStockData(symbol));

    // Daily time step
    double dt = 1.0 / 252;

    std::vector<std::vector<double>> simulations(numSimulations);
    for (int i = 0; i < numSimulations; ++i) {
        simulations[i] = geometricBrownianMotion(currentPrice, mu, sigma, timeSteps, dt);
    }

    return simulations;
}