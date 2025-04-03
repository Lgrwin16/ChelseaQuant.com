#ifndef IB_CONNECTOR_H
#define IB_CONNECTOR_H

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <queue>

// Interactive Brokers API headers
#include <EClient.h>
#include <EWrapper.h>
#include <EReaderOSSignal.h>
#include <EReader.h>
#include <Contract.h>
#include <Order.h>
#include <Bar.h>

#include "MarketData.h"

struct TradeOrder {
    std::string symbol;
    int quantity;
    double price;
    bool isBuy;
    std::string orderType; // "MARKET", "LIMIT", etc.
    std::string timeInForce; // "DAY", "GTC", etc.
};

// Forward declaration
class EReaderGuard;

class IBConnector : public EWrapper, public EClient {
private:
    // IB API components
    std::unique_ptr<EReaderOSSignal> m_osSignal;
    std::unique_ptr<EReader> m_reader;
    std::atomic<bool> m_isConnected;
    std::atomic<int> m_nextOrderId;
    std::atomic<int> m_nextRequestId;

    // Connection parameters
    std::string m_host;
    int m_port;
    int m_clientId;

    // Synchronization
    std::mutex m_mutex;
    std::condition_variable m_cv;
    bool m_waitingForData;

    // Data storage
    std::map<int, PriceTimeSeries> m_historicalData;
    std::map<std::string, PriceData> m_latestPrices;
    std::map<int, TradeOrder> m_openOrders;
    std::map<std::string, double> m_positions;
    double m_accountCash;

    // Request to symbol mapping
    std::map<int, std::string> m_requestToSymbol;

    // Order tracking
    std::map<int, bool> m_orderPlaced;

    // Helper methods
    int getNextRequestId();
    Contract createContract(const std::string& symbol);
    Order createOrder(const TradeOrder& orderDetails);

public:
    IBConnector(const std::string& host = "127.0.0.1",
               int port = 7496,
               int clientId = 1);
    ~IBConnector();

    bool connect();
    void disconnect();
    bool isConnected() const;
    void processMessages();

    // Market data methods
    bool requestMarketData(const std::vector<std::string>& symbols);
    PriceData getLatestPrice(const std::string& symbol) const;
    PriceTimeSeries getHistoricalData(const std::string& symbol,
                                     const std::string& endDateTime,
                                     int duration,
                                     const std::string& barSize,
                                     const std::string& whatToShow);

    // Trading methods
    int placeOrder(const TradeOrder& order);
    bool cancelOrder(int orderId);
    std::map<int, TradeOrder> getOpenOrders() const;

    // Account methods
    double getAccountCash() const;
    double getPositionSize(const std::string& symbol) const;
    std::vector<std::string> getPositions() const;
    void requestAccountUpdates();

    // EWrapper interface implementation
    void tickPrice(TickerId tickerId, TickType field, double price, const TickAttrib& attrib) override;
    void tickSize(TickerId tickerId, TickType field, int size) override;
    void orderStatus(OrderId orderId, const std::string& status, double filled,
                    double remaining, double avgFillPrice, int permId, int parentId,
                    double lastFillPrice, int clientId, const std::string& whyHeld, double mktCapPrice) override;
    void openOrder(OrderId orderId, const Contract&, const Order&, const OrderState&) override;
    void historicalData(TickerId reqId, const Bar& bar) override;
    void historicalDataEnd(int reqId, const std::string& startDateStr, const std::string& endDateStr) override;
    void error(int id, int errorCode, const std::string& errorString) override;
    void connectionClosed() override;
    void nextValidId(OrderId orderId) override;
    void position(const std::string& account, const Contract& contract, double position, double avgCost) override;
    void positionEnd() override;
    void accountSummary(int reqId, const std::string& account, const std::string& tag, const std::string& value, const std::string& currency) override;
    void accountSummaryEnd(int reqId) override;
};

#endif // IB_CONNECTOR_H