cmake_minimum_required(VERSION 3.12)

project(matchtest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(flann REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(NLopt REQUIRED)

add_executable(matchLogic src/main.cpp
						include/line_fit.h
						src/line_fit.cpp
						src/method_map.cpp 
						src/method_flann.cpp 
						src/method_fitCurve.cpp
						src/method_getDistance.cpp)
target_include_directories(matchLogic PRIVATE include)
target_link_libraries(matchLogic Eigen3::Eigen flann::flann NLopt::nlopt)
