#include "IBConnector.h"
#include <iostream>
#include <thread>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <ctime>

IBConnector::IBConnector(const std::string& host, int port, int clientId)
    : EClient(this, nullptr),
      m_host(host),
      m_port(port),
      m_clientId(clientId),
      m_isConnected(false),
      m_nextOrderId(-1),
      m_nextRequestId(1),
      m_waitingForData(false),
      m_accountCash(0.0) {

    m_osSignal = std::make_unique<EReaderOSSignal>(2000); // 2-second timeout
}

IBConnector::~IBConnector() {
    disconnect();
}

bool IBConnector::connect() {
    if (m_isConnected) {
        std::cout << "Already connected to Interactive Brokers" << std::endl;
        return true;
    }

    std::cout << "Connecting to Interactive Brokers at " << m_host << ":" << m_port
              << " with client ID " << m_clientId << std::endl;

    // Connect to TWS or IB Gateway
    EClient::connect(m_host.c_str(), m_port, m_clientId);

    // Wait for nextValidId callback, which confirms connection
    std::unique_lock<std::mutex> lock(m_mutex);
    m_waitingForData = true;
    if (!m_cv.wait_for(lock, std::chrono::seconds(10), [this] { return m_nextOrderId >= 0; })) {
        std::cerr << "Connection to Interactive Brokers timed out" << std::endl;
        EClient::disconnect();
        m_waitingForData = false;
        return false;
    }

    m_waitingForData = false;
    m_isConnected = true;

    // Start message processing in a separate thread
    m_reader = std::make_unique<EReader>(this, m_osSignal.get());
    m_reader->start();

    std::thread([this]() {
        while (m_isConnected) {
            processMessages();
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }).detach();

    std::cout << "Connected to Interactive Brokers" << std::endl;
    return true;
}

void IBConnector::disconnect() {
    if (!m_isConnected) {
        return;
    }

    std::cout << "Disconnecting from Interactive Brokers" << std::endl;

    EClient::disconnect();
    m_isConnected = false;
}

bool IBConnector::isConnected() const {
    return m_isConnected && EClient::isConnected();
}

void IBConnector::processMessages() {
    if (!m_isConnected) return;

    m_osSignal->waitForSignal();
    m_reader->processMsgs();
}

bool IBConnector::requestMarketData(const std::vector<std::string>& symbols) {
    if (!isConnected()) {
        std::cerr << "Not connected to Interactive Brokers" << std::endl;
        return false;
    }

    std::cout << "Requesting market data for " << symbols.size() << " symbols" << std::endl;

    for (const auto& symbol : symbols) {
        Contract contract = createContract(symbol);
        int reqId = getNextRequestId();

        m_requestToSymbol[reqId] = symbol;
        EClient::reqMktData(reqId, contract, "", false, false, TagValueListSPtr());

        std::cout << "  Subscribing to market data for " << symbol << " (reqId: " << reqId << ")" << std::endl;
    }

    return true;
}

PriceData IBConnector::getLatestPrice(const std::string& symbol) const {
    std::lock_guard<std::mutex> lock(const_cast<IBConnector*>(this)->m_mutex);

    auto it = m_latestPrices.find(symbol);
    if (it == m_latestPrices.end()) {
        throw std::runtime_error("No price data available for symbol: " + symbol);
    }

    return it->second;
}

PriceTimeSeries IBConnector::getHistoricalData(
    const std::string& symbol,
    const std::string& endDateTime,
    int duration,
    const std::string& barSize,
    const std::string& whatToShow) {

    if (!isConnected()) {
        throw std::runtime_error("Not connected to Interactive Brokers");
    }

    std::cout << "Fetching historical data for " << symbol
              << " ending at " << endDateTime
              << " with duration " << duration
              << " and bar size " << barSize << std::endl;

    // Create contract
    Contract contract = createContract(symbol);

    // Format duration string (e.g., "30 D" for 30 days)
    std::string durationStr = std::to_string(duration) + " D";

    // Request historical data
    int reqId = getNextRequestId();
    m_requestToSymbol[reqId] = symbol;

    // Clear any previous data for this request
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_historicalData[reqId].clear();
        m_waitingForData = true;
    }

    EClient::reqHistoricalData(
        reqId,
        contract,
        endDateTime,
        durationStr,
        barSize,
        whatToShow,
        1, // useRTH = 1 for regular trading hours only
        1, // formatDate = 1 for yyyyMMdd format
        false,
        TagValueListSPtr()
    );

    // Wait for data to be received
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        if (!m_cv.wait_for(lock, std::chrono::seconds(30), [this, reqId] {
            return !m_waitingForData || !m_historicalData[reqId].empty();
        })) {
            throw std::runtime_error("Timeout waiting for historical data for " + symbol);
        }
    }

    // Return the collected data
    std::lock_guard<std::mutex> lock(m_mutex);
    return m_historicalData[reqId];
}

int IBConnector::placeOrder(const TradeOrder& orderDetails) {
    if (!isConnected()) {
        throw std::runtime_error("Not connected to Interactive Brokers");
    }

    if (m_nextOrderId < 0) {
        throw std::runtime_error("No valid order ID available");
    }

    // Create contract
    Contract contract = createContract(orderDetails.symbol);

    // Create order
    Order order = createOrder(orderDetails);

    // Get next order ID
    int orderId = m_nextOrderId++;

    // Place order
    std::cout << "Placing " << (orderDetails.isBuy ? "BUY" : "SELL") << " order for " << orderDetails.symbol
              << ", Quantity: " << orderDetails.quantity
              << ", Price: " << orderDetails.price
              << ", Type: " << orderDetails.orderType
              << ", OrderId: " << orderId << std::endl;

    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_orderPlaced[orderId] = false;
        m_openOrders[orderId] = orderDetails;
    }

    EClient::placeOrder(orderId, contract, order);

    // Wait for order confirmation
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        if (!m_cv.wait_for(lock, std::chrono::seconds(5), [this, orderId] {
            return m_orderPlaced[orderId];
        })) {
            std::cerr << "Warning: No immediate confirmation for order " << orderId << std::endl;
        }
    }

    return orderId;
}

bool IBConnector::cancelOrder(int orderId) {
    if (!isConnected()) {
        throw std::runtime_error("Not connected to Interactive Brokers");
    }

    std::cout << "Cancelling order " << orderId << std::endl;

    EClient::cancelOrder(orderId, "");

    // Remove from open orders (actual removal will happen on callback)
    return true;
}

std::map<int, TradeOrder> IBConnector::getOpenOrders() const {
    std::lock_guard<std::mutex> lock(const_cast<IBConnector*>(this)->m_mutex);
    return m_openOrders;
}

double IBConnector::getAccountCash() const {
    std::lock_guard<std::mutex> lock(const_cast<IBConnector*>(this)->m_mutex);
    return m_accountCash;
}

double IBConnector::getPositionSize(const std::string& symbol) const {
    std::lock_guard<std::mutex> lock(const_cast<IBConnector*>(this)->m_mutex);

    auto it = m_positions.find(symbol);
    if (it == m_positions.end()) {
        return 0.0;  // No position
    }

    return it->second;
}

std::vector<std::string> IBConnector::getPositions() const {
    std::lock_guard<std::mutex> lock(const_cast<IBConnector*>(this)->m_mutex);

    std::vector<std::string> symbolsList;
    for (const auto& position : m_positions) {
        if (position.second != 0.0) {
            symbolsList.push_back(position.first);
        }
    }

    return symbolsList;
}

void IBConnector::requestAccountUpdates() {
    if (!isConnected()) {
        throw std::runtime_error("Not connected to Interactive Brokers");
    }

    // Request account summary
    int reqId = getNextRequestId();
    EClient::reqAccountSummary(reqId, "All", "TotalCashValue");

    // Request positions
    EClient::reqPositions();
}

int IBConnector::getNextRequestId() {
    return m_nextRequestId++;
}

Contract IBConnector::createContract(const std::string& symbol) {
    Contract contract;
    contract.symbol = symbol;
    contract.secType = "STK";
    contract.currency = "USD";
    contract.exchange = "SMART";

    return contract;
}

Order IBConnector::createOrder(const TradeOrder& orderDetails) {
    Order order;

    // Set action (BUY/SELL)
    order.action = orderDetails.isBuy ? "BUY" : "SELL";

    // Set quantity
    order.totalQuantity = orderDetails.quantity;

    // Set order type
    if (orderDetails.orderType == "MARKET") {
        order.orderType = "MKT";
    } else if (orderDetails.orderType == "LIMIT") {
        order.orderType = "LMT";
        order.lmtPrice = orderDetails.price;
    } else {
        order.orderType = orderDetails.orderType;
    }

    // Set time in force
    if (orderDetails.timeInForce == "DAY") {
        order.tif = "DAY";
    } else if (orderDetails.timeInForce == "GTC") {
        order.tif = "GTC";
    } else {
        order.tif = orderDetails.timeInForce;
    }

    // Additional order properties
    order.outsideRth = false;  // Regular trading hours only
    order.transmit = true;     // Transmit order immediately

    return order;
}

// EWrapper interface implementations

void IBConnector::tickPrice(TickerId tickerId, TickType field, double price, const TickAttrib& attrib) {
    std::lock_guard<std::mutex> lock(m_mutex);

    auto it = m_requestToSymbol.find(tickerId);
    if (it == m_requestToSymbol.end()) {
        return;
    }

    std::string symbol = it->second;

    // Update price data based on tick type
    if (field == TickType::LAST || field == TickType::DELAYED_LAST) {
        m_latestPrices[symbol].close = price;
    } else if (field == TickType::BID || field == TickType::DELAYED_BID) {
        // Not used in this example but could be stored
    } else if (field == TickType::ASK || field == TickType::DELAYED_ASK) {
        // Not used in this example but could be stored
    } else if (field == TickType::HIGH || field == TickType::DELAYED_HIGH) {
        m_latestPrices[symbol].high = price;
    } else if (field == TickType::LOW || field == TickType::DELAYED_LOW) {
        m_latestPrices[symbol].low = price;
    } else if (field == TickType::OPEN || field == TickType::DELAYED_OPEN) {
        m_latestPrices[symbol].open = price;
    }

    // Add timestamp
    auto now = std::time(nullptr);
    auto tm = *std::localtime(&now);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    m_latestPrices[symbol].timestamp = oss.str();
}

void IBConnector::tickSize(TickerId tickerId, TickType field, int size) {
    std::lock_guard<std::mutex> lock(m_mutex);

    auto it = m_requestToSymbol.find(tickerId);
    if (it == m_requestToSymbol.end()) {
        return;
    }

    std::string symbol = it->second;

    if (field == TickType::VOLUME || field == TickType::DELAYED_VOLUME) {
        m_latestPrices[symbol].volume = size;
    }
}

void IBConnector::orderStatus(OrderId orderId, const std::string& status, double filled,
                            double remaining, double avgFillPrice, int permId, int parentId,
                            double lastFillPrice, int clientId, const std::string& whyHeld, double mktCapPrice) {
    std::lock_guard<std::mutex> lock(m_mutex);

    std::cout << "Order status update - ID: " << orderId
              << ", Status: " << status
              << ", Filled: " << filled
              << ", Remaining: " << remaining
              << ", Avg Fill Price: " << avgFillPrice << std::endl;

    if (status == "Filled") {
        // Remove from open orders if completely filled
        if (remaining == 0) {
            m_openOrders.erase(orderId);
        }
    } else if (status == "Cancelled" || status == "Inactive") {
        m_openOrders.erase(orderId);
    }

    // Mark order as placed (for wait condition)
    m_orderPlaced[orderId] = true;
    m_cv.notify_all();
}

void IBConnector::openOrder(OrderId orderId, const Contract& contract, const Order& order, const OrderState& orderState) {
    std::lock_guard<std::mutex> lock(m_mutex);

    std::cout << "Open order notification - ID: " << orderId
              << ", Symbol: " << contract.symbol
              << ", Action: " << order.action
              << ", Quantity: " << order.totalQuantity
              << ", Status: " << orderState.status << std::endl;

    // Update or add to open orders map
    TradeOrder tradeOrder;
    tradeOrder.symbol = contract.symbol;
    tradeOrder.quantity = order.totalQuantity;
    tradeOrder.isBuy = (order.action == "BUY");
    tradeOrder.orderType = order.orderType;
    tradeOrder.timeInForce = order.tif;

    if (order.orderType == "LMT") {
        tradeOrder.price = order.lmtPrice;
    }

    m_openOrders[orderId] = tradeOrder;
}

void IBConnector::historicalData(TickerId reqId, const Bar& bar) {
    std::lock_guard<std::mutex> lock(m_mutex);

    auto it = m_requestToSymbol.find(reqId);
    if (it == m_requestToSymbol.end()) {
        return;
    }

    PriceData data;
    data.open = bar.open;
    data.high = bar.high;
    data.low = bar.low;
    data.close = bar.close;
    data.adjustedClose = bar.close;  // No adjusted close in IB data
    data.volume = bar.volume;
    data.timestamp = bar.time;

    m_historicalData[reqId].push_back(data);
}

void IBConnector::historicalDataEnd(int reqId, const std::string& startDateStr, const std::string& endDateStr) {
    std::lock_guard<std::mutex> lock(m_mutex);
    std::cout << "Historical data received for reqId: " << reqId
              << ", from " << startDateStr
              << " to " << endDateStr
              << ", " << m_historicalData[reqId].size() << " bars" << std::endl;

    m_waitingForData = false;
    m_cv.notify_all();
}

void IBConnector::error(int id, int errorCode, const std::string& errorString) {
    std::cerr << "IB Error " << errorCode << ": " << errorString;

    if (id > 0) {
        std::cerr << " (ID: " << id << ")";
    }

    std::cerr << std::endl;

    // Handle specific error codes
    if (errorCode == 502) {  // Cannot connect to TWS
        m_isConnected = false;
    } else if (id == -1 && errorCode == 2104) {  // Market data farm connection is OK
        // Informational message, not an error
    } else if (id == -1 && errorCode == 2106) {  // HMDS data farm connection is OK
        // Informational message, not an error
    }

    if (m_waitingForData) {
        m_waitingForData = false;
        m_cv.notify_all();
    }
}

void IBConnector::connectionClosed() {
    std::cout << "Connection to Interactive Brokers closed" << std::endl;
    m_isConnected = false;

    if (m_waitingForData) {
        m_waitingForData = false;
        m_cv.notify_all();
    }
}

void IBConnector::nextValidId(OrderId orderId) {
    std::lock_guard<std::mutex> lock(m_mutex);

    std::cout << "Next valid order ID: " << orderId << std::endl;
    m_nextOrderId = orderId;

    m_cv.notify_all();  // Signal connection is ready
}

void IBConnector::position(const std::string& account, const Contract& contract, double position, double avgCost) {
    std::lock_guard<std::mutex> lock(m_mutex);

    std::cout << "Position update - Account: " << account
              << ", Symbol: " << contract.symbol
              << ", Position: " << position
              << ", Avg Cost: " << avgCost << std::endl;

    m_positions[contract.symbol] = position;
}

void IBConnector::positionEnd() {
    std::cout << "Position updates completed" << std::endl;
}

void IBConnector::accountSummary(int reqId, const std::string& account, const std::string& tag, const std::string& value, const std::string& currency) {
    std::lock_guard<std::mutex> lock(m_mutex);

    std::cout << "Account summary - Account: " << account
              << ", Tag: " << tag
              << ", Value: " << value
              << ", Currency: " << currency << std::endl;

    if (tag == "TotalCashValue") {
        try {
            m_accountCash = std::stod(value);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing cash value: " << e.what() << std::endl;
        }
    }
}

void IBConnector::accountSummaryEnd(int reqId) {
    std::cout << "Account summary updates completed" << std::endl;
}