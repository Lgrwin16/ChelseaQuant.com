#ifndef MARKET_DATA_H
#define MARKET_DATA_H

#include <string>
#include <vector>
#include <map>

struct PriceData {
    double open;
    double high;
    double low;
    double close;
    double adjustedClose;
    long volume;
    std::string timestamp;
};

typedef std::vector<PriceData> PriceTimeSeries;

class MarketData {
private:
    std::map<std::string, PriceTimeSeries> stockData;
    PriceTimeSeries marketIndexData;
    std::string marketIndex;
    int lookbackDays;

public:
    MarketData(const std::string& marketIndex, int lookbackDays);
    ~MarketData();

    void fetchHistoricalData(const std::vector<std::string>& symbols);
    void updateLiveData(const std::vector<std::string>& symbols);

    PriceTimeSeries getStockData(const std::string& symbol) const;
    PriceTimeSeries getMarketIndexData() const;

    std::vector<double> getReturns(const std::string& symbol) const;
    std::vector<double> getMarketReturns() const;

    double getCurrentPrice(const std::string& symbol) const;
    double getMovingAverage(const std::string& symbol, int period) const;
    double getVolatility(const std::string& symbol, int period) const;
};

#endif // MARKET_DATA_H