
cmake_minimum_required(VERSION 3.1)

project(poissonCos)

# Set project sources
set(SRCS
     poissonCos.cpp
     prediction-projection/Prediction.cpp
     prediction-projection/Projection.cpp
     )
find_package(NEOS)
add_executable(poissonCos ${SRCS})
target_link_libraries(poissonCos neos::neos)

