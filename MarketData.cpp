#include "MarketData.h"
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <cmath>

MarketData::MarketData(const std::string& marketIndex, int lookbackDays)
    : marketIndex(marketIndex), lookbackDays(lookbackDays) {
}

MarketData::~MarketData() {
    // Clean up resources if needed
}

void MarketData::fetchHistoricalData(const std::vector<std::string>& symbols) {
    // In a real implementation, this would use the IB API to fetch historical data
    // For demonstration purposes, we'll simulate it

    std::cout << "Fetching historical data for " << symbols.size() << " symbols" << std::endl;

    // Clear existing data
    stockData.clear();
    marketIndexData.clear();

    // For each symbol, generate placeholder data
    for (const auto& symbol : symbols) {
        PriceTimeSeries timeSeries;

        // Generate lookbackDays days of data
        for (int i = 0; i < lookbackDays; ++i) {
            PriceData data;

            // In a real implementation, these would be actual market data
            double basePrice = (symbol == "AAPL") ? 180.0 :
                              (symbol == "MSFT") ? 350.0 :
                              (symbol == "GOOG") ? 140.0 :
                              (symbol == "SPY") ? 450.0 : 100.0;

            // Generate some random price movement
            double randomFactor = 0.98 + (std::rand() % 5) / 100.0;
            double dayPrice = basePrice * randomFactor * (1.0 + 0.001 * (lookbackDays - i));

            data.open = dayPrice * 0.995;
            data.high = dayPrice * 1.01;
            data.low = dayPrice * 0.99;
            data.close = dayPrice;
            data.adjustedClose = dayPrice;
            data.volume = 1000000 + std::rand() % 1000000;

            // Format: YYYY-MM-DD
            char dateBuffer[11];
            std::sprintf(dateBuffer, "2023-%02d-%02d", 1 + i / 30, 1 + (i % 30));
            data.timestamp = dateBuffer;

            timeSeries.push_back(data);
        }

        if (symbol == marketIndex) {
            marketIndexData = timeSeries;
        } else {
            stockData[symbol] = timeSeries;
        }
    }

    std::cout << "Historical data fetched successfully" << std::endl;
}

void MarketData::updateLiveData(const std::vector<std::string>& symbols) {
    // In a real implementation, this would use IB API to fetch live data
    // For demonstration purposes, we'll simulate live data

    std::cout << "Updating live data for " << symbols.size() << " symbols" << std::endl;

    for (const auto& symbol : symbols) {
        // Get the last historical data point
        PriceTimeSeries& series = (symbol == marketIndex) ? marketIndexData : stockData[symbol];

        if (series.empty()) {
            std::cerr << "No historical data available for " << symbol << std::endl;
            continue;
        }

        PriceData lastData = series.back();

        // Simulate a new data point with slight randomness
        PriceData newData;
        double randomFactor = 0.995 + (std::rand() % 11) / 1000.0;  // Random between -0.5% and 0.5%

        newData.close = lastData.close * randomFactor;
        newData.open = newData.close * 0.998;
        newData.high = newData.close * 1.005;
        newData.low = newData.close * 0.995;
        newData.adjustedClose = newData.close;
        newData.volume = lastData.volume * (0.9 + (std::rand() % 21) / 100.0);  // Random volume ±10%

        // Current timestamp
        newData.timestamp = "LIVE";

        // Add to beginning of series (most recent first)
        series.insert(series.begin(), newData);

        // Keep only lookbackDays + 1 data points (including live data)
        if (series.size() > lookbackDays + 1) {
            series.resize(lookbackDays + 1);
        }
    }

    std::cout << "Live data updated successfully" << std::endl;
}

PriceTimeSeries MarketData::getStockData(const std::string& symbol) const {
    auto it = stockData.find(symbol);
    if (it == stockData.end()) {
        throw std::runtime_error("No data available for symbol: " + symbol);
    }
    return it->second;
}

PriceTimeSeries MarketData::getMarketIndexData() const {
    return marketIndexData;
}

std::vector<double> MarketData::getReturns(const std::string& symbol) const {
    std::vector<double> returns;

    PriceTimeSeries series;
    try {
        series = getStockData(symbol);
    } catch (const std::exception& e) {
        if (symbol == marketIndex) {
            series = getMarketIndexData();
        } else {
            throw;
        }
    }

    // Calculate daily returns
    for (size_t i = 1; i < series.size(); ++i) {
        double prevClose = series[i].adjustedClose;
        double currentClose = series[i-1].adjustedClose;  // i-1 because recent data is at the beginning

        // Calculate log return
        double logReturn = std::log(currentClose / prevClose);
        returns.push_back(logReturn);
    }

    return returns;
}

std::vector<double> MarketData::getMarketReturns() const {
    return getReturns(marketIndex);
}

double MarketData::getCurrentPrice(const std::string& symbol) const {
    PriceTimeSeries series;
    try {
        series = getStockData(symbol);
    } catch (const std::exception& e) {
        if (symbol == marketIndex) {
            series = getMarketIndexData();
        } else {
            throw;
        }
    }

    if (series.empty()) {
        throw std::runtime_error("No price data available for symbol: " + symbol);
    }

    // Return the most recent close price (at index 0)
    return series[0].close;
}

double MarketData::getMovingAverage(const std::string& symbol, int period) const {
    PriceTimeSeries series;
    try {
        series = getStockData(symbol);
    } catch (const std::exception& e) {
        if (symbol == marketIndex) {
            series = getMarketIndexData();
        } else {
            throw;
        }
    }

    if (series.size() < period) {
        throw std::runtime_error("Insufficient data for calculating moving average");
    }

    double sum = 0.0;
    for (int i = 0; i < period; ++i) {
        sum += series[i].close;
    }

    return sum / period;
}

double MarketData::getVolatility(const std::string& symbol, int period) const {
    std::vector<double> returns = getReturns(symbol);

    if (returns.size() < period) {
        throw std::runtime_error("Insufficient data for calculating volatility");
    }

    // Calculate mean of returns
    double sum = 0.0;
    for (int i = 0; i < period; ++i) {
        sum += returns[i];
    }
    double mean = sum / period;

    // Calculate variance
    double sumSquaredDiff = 0.0;
    for (int i = 0; i < period; ++i) {
        double diff = returns[i] - mean;
        sumSquaredDiff += diff * diff;
    }
    double variance = sumSquaredDiff / period;

    // Return annualized volatility (assuming daily returns)
    return std::sqrt(variance) * std::sqrt(252);  // 252 trading days in a year
}