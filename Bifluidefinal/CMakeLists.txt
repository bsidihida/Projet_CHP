
cmake_minimum_required(VERSION 3.1)

project(TaylorGreen)

# Set project sources
set(SRCS
   #  prediction-projection/Gradient2.cpp
     TaylorGreen.cpp
     prediction-projection/Prediction.cpp
     prediction-projection/Projection.cpp
    
     )
find_package(NEOS)
add_executable(TaylorGreen ${SRCS})
target_link_libraries(TaylorGreen neos::neos)

