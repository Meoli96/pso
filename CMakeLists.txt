cmake_minimum_required(VERSION 3.30)
project(pso)

set (CMAKE_CXX_STANDARD 20)

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

add_executable(pso ${SRC_DIR}/main.cpp)
# Add another executable compiled with dod flag


add_executable(pso_dod ${SRC_DIR}/main.cpp)
target_compile_definitions(pso_dod PRIVATE DOD_TEST)

