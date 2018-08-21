cmake_minimum_required(VERSION 3.6)
project(train_schedule_optimisation)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(Cplex)
include_directories(SYSTEM ${CPLEX_INCLUDE_DIRS})

set(SOURCE_FILES main.cpp Instance.h Train.cpp Train.h Requirement.h Route.h)
add_executable(train_schedule_optimisation ${SOURCE_FILES})
target_link_libraries(train_schedule_optimisation ${CPLEX_LIBRARIES})
