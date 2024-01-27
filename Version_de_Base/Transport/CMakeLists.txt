
cmake_minimum_required(VERSION 3.1)

project(transportAtmospheric)

# Set project sources
set(SRCS
     transportAtmospheric.cpp
     prediction-projection/Prediction.cpp
     prediction-projection/Projection.cpp
     )
find_package(NEOS)
add_executable(transportAtmospheric ${SRCS})
target_link_libraries(transportAtmospheric neos::neos)

