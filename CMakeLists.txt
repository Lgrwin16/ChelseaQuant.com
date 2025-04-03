cmake_minimum_required(VERSION 3.10)
project(AlgorithmicTradingSystem)

set(CMAKE_CXX_STANDARD 17)

# Find IBAPI location - customize this path based on your setup
set(IBAPI_DIR "${CMAKE_SOURCE_DIR}/ibapi")

# Check if IBAPI directory exists
if(NOT EXISTS ${IBAPI_DIR})
    message(FATAL_ERROR "Interactive Brokers API directory not found at ${IBAPI_DIR}. Please update the IBAPI_DIR variable in CMakeLists.txt.")
endif()

# IBAPI include directories
include_directories(
        ${IBAPI_DIR}/client
        ${IBAPI_DIR}/source/cppclient/client
        ${IBAPI_DIR}/source/cppclient/posix
        ${CMAKE_SOURCE_DIR}    # For local header files
)

# Explicitly list all source files
set(SOURCE_FILES
        # Main application file
        main.cpp

        # Core algorithm components
        TradingAlgorithm.cpp
        MarketData.cpp
        Statistics.cpp
        IBConnector.cpp
        TradeExecution.cpp
)

# Get list of IB API source files
file(GLOB IBAPI_SOURCE_FILES
        "${IBAPI_DIR}/source/cppclient/client/*.cpp"
        "${IBAPI_DIR}/source/cppclient/posix/*.cpp"
)

# Add executable with all source files
add_executable(AlgorithmicTradingSystem
        ${SOURCE_FILES}
        ${IBAPI_SOURCE_FILES}
)

# Link with pthread on Unix-like systems
if(UNIX)
    target_link_libraries(AlgorithmicTradingSystem pthread)
endif()

# Set additional compiler flags if needed
if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
else()
    add_compile_options(-Wall -Wextra)
endif()

# Install target
install(TARGETS AlgorithmicTradingSystem
        RUNTIME DESTINATION bin)

# Enable testing
enable_testing()