#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include "TradingAlgorithm.h"
#include "TradeExecution.h"

void printUsage() {
    std::cout << "Usage: TradingAlgorithm <strategy> <symbols> <marketIndex> <riskFreeRate> [options]" << std::endl;
    std::cout << "  strategy     - STAT_ARB or MEAN_REVERSION" << std::endl;
    std::cout << "  symbols      - Comma-separated list of stock symbols" << std::endl;
    std::cout << "  marketIndex  - Market index symbol (e.g., SPY)" << std::endl;
    std::cout << "  riskFreeRate - Risk-free rate (e.g., 0.03 for 3%)" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  --lookback N - Lookback period in days (default: 252)" << std::endl;
    std::cout << "  --threshold N - Confidence threshold (default: 1.5)" << std::endl;
    std::cout << "  --backtest START END - Run backtest from START to END dates (YYYY-MM-DD)" << std::endl;
    std::cout << "  --position-size N - Maximum position size as percentage (default: 5%)" << std::endl;
    std::cout << "  --risk-per-trade N - Risk per trade as percentage (default: 1%)" << std::endl;
    std::cout << "  --position-model MODEL - Position sizing model (FIXED, VOLATILITY, KELLY, EQUAL_RISK, VAR)" << std::endl;
    std::cout << "  --initial-capital N - Initial capital in USD (default: 1000000)" << std::endl;
}

std::vector<std::string> parseSymbols(const std::string& symbolsStr) {
    std::vector<std::string> symbols;
    std::string symbol;
    size_t pos = 0;
    std::string str = symbolsStr;

    while ((pos = str.find(',')) != std::string::npos) {
        symbol = str.substr(0, pos);
        symbols.push_back(symbol);
        str.erase(0, pos + 1);
    }
    symbols.push_back(str);  // Add the last symbol

    return symbols;
}

PositionSizingModel parsePositionModel(const std::string& modelStr) {
    if (modelStr == "FIXED") return PositionSizingModel::FIXED_PERCENTAGE;
    if (modelStr == "VOLATILITY") return PositionSizingModel::VOLATILITY_ADJUSTED;
    if (modelStr == "KELLY") return PositionSizingModel::KELLY_CRITERION;
    if (modelStr == "EQUAL_RISK") return PositionSizingModel::EQUAL_RISK;
    if (modelStr == "VAR") return PositionSizingModel::VALUE_AT_RISK;

    // Default
    return PositionSizingModel::VOLATILITY_ADJUSTED;
}

int main(int argc, char* argv[]) {
    // Default parameters
    TradingStrategy strategy = TradingStrategy::STAT_ARB;
    std::vector<std::string> symbols = {"AAPL", "MSFT", "GOOG", "AMZN", "META"};
    std::string marketIndex = "SPY";
    double riskFreeRate = 0.03;
    double lookbackPeriod = 252;
    double confidenceThreshold = 1.5;
    double initialCapital = 1000000.0;

    // Position sizing defaults
    PositionSizing positionSizing;
    positionSizing.model = PositionSizingModel::VOLATILITY_ADJUSTED;
    positionSizing.maxPositionSize = 0.05;  // 5% max position size
    positionSizing.riskPerTrade = 0.01;     // 1% risk per trade

    // Parse command line arguments
    if (argc >= 5) {
        // Parse required arguments
        std::string strategyStr = argv[1];
        strategy = (strategyStr == "MEAN_REVERSION") ? TradingStrategy::MEAN_REVERSION : TradingStrategy::STAT_ARB;

        symbols = parseSymbols(argv[2]);
        marketIndex = argv[3];
        riskFreeRate = std::stod(argv[4]);

        // Parse optional arguments
        for (int i = 5; i < argc; i++) {
            std::string arg = argv[i];

            if (arg == "--lookback" && i + 1 < argc) {
                lookbackPeriod = std::stod(argv[++i]);
            }
            else if (arg == "--threshold" && i + 1 < argc) {
                confidenceThreshold = std::stod(argv[++i]);
            }
            else if (arg == "--position-size" && i + 1 < argc) {
                positionSizing.maxPositionSize = std::stod(argv[++i]) / 100.0;
            }
            else if (arg == "--risk-per-trade" && i + 1 < argc) {
                positionSizing.riskPerTrade = std::stod(argv[++i]) / 100.0;
            }
            else if (arg == "--position-model" && i + 1 < argc) {
                positionSizing.model = parsePositionModel(argv[++i]);
            }
            else if (arg == "--initial-capital" && i + 1 < argc) {
                initialCapital = std::stod(argv[++i]);
            }
            else if (arg == "--backtest" && i + 2 < argc) {
                std::string startDate = argv[++i];
                std::string endDate = argv[++i];
                // Backtest handling would be added here
            }
        }
    } else {
        printUsage();
        return 1;
    }

    // Print configuration
    std::cout << "Trading Algorithm Configuration:" << std::endl;
    std::cout << "  Strategy: " << (strategy == TradingStrategy::STAT_ARB ? "Statistical Arbitrage" : "Mean Reversion") << std::endl;
    std::cout << "  Symbols: ";
    for (const auto& symbol : symbols) {
        std::cout << symbol << " ";
    }
    std::cout << std::endl;
    std::cout << "  Market Index: " << marketIndex << std::endl;
    std::cout << "  Risk-Free Rate: " << riskFreeRate * 100 << "%" << std::endl;
    std::cout << "  Lookback Period: " << lookbackPeriod << " days" << std::endl;
    std::cout << "  Confidence Threshold: " << confidenceThreshold << std::endl;
    std::cout << "  Initial Capital: $" << initialCapital << std::endl;
    std::cout << "  Position Sizing: ";
    switch (positionSizing.model) {
        case PositionSizingModel::FIXED_PERCENTAGE:
            std::cout << "Fixed Percentage";
            break;
        case PositionSizingModel::VOLATILITY_ADJUSTED:
            std::cout << "Volatility Adjusted";
            break;
        case PositionSizingModel::KELLY_CRITERION:
            std::cout << "Kelly Criterion";
            break;
        case PositionSizingModel::EQUAL_RISK:
            std::cout << "Equal Risk";
            break;
        case PositionSizingModel::VALUE_AT_RISK:
            std::cout << "Value at Risk";
            break;
    }
    std::cout << " (Max: " << positionSizing.maxPositionSize * 100 << "%, Risk: " << positionSizing.riskPerTrade * 100 << "%)" << std::endl;

    try {
        // Create trading algorithm
        TradingAlgorithm algorithm(strategy, symbols, marketIndex, riskFreeRate, lookbackPeriod, confidenceThreshold);

        // Initialize algorithm
        algorithm.initialize();

        // Create IB connector reference
        std::shared_ptr<IBConnector> ibConnector = std::make_shared<IBConnector>();
        if (!ibConnector->connect()) {
            std::cerr << "Failed to connect to Interactive Brokers" << std::endl;
            return 1;
        }

        // Create market data reference
        std::shared_ptr<MarketData> marketData = std::make_shared<MarketData>(marketIndex, static_cast<int>(lookbackPeriod));

        // Create statistics reference
        std::shared_ptr<Statistics> statistics = std::make_shared<Statistics>(marketData);

        // Create trade execution module
        TradeExecution tradeExecution(ibConnector, statistics, initialCapital, positionSizing);

        // Fetch and analyze data
        algorithm.fetchData();
        algorithm.analyzeData();

        // Generate trading signals
        algorithm.generateTrades();

        // Execute trades (this would be handled through the trade execution module in a real system)
        algorithm.executeTrades();

        // Report execution status
        tradeExecution.updateAccountInfo();
        std::cout << "Portfolio Value: $" << tradeExecution.getTotalPortfolioValue() << std::endl;
        std::cout << "Portfolio Heat: " << tradeExecution.calculatePortfolioHeat() * 100 << "%" << std::endl;

        // Clean up
        ibConnector->disconnect();
        std::cout << "Trading algorithm completed successfully" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}