cmake_minimum_required(VERSION 3.6)
project(train_schedule_optimisation)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp Instance.h Train.cpp Train.h Requirement.h Route.h)
add_executable(train_schedule_optimisation ${SOURCE_FILES})