cmake_minimum_required(VERSION 3.0)
project(dork_examples)

set(CMAKE_CXX_STANDARD 17)
set(PROJECT_DIR ${PROJECT_SOURCE_DIR}/../../)

include_directories(${PROJECT_DIR}/include)

file(GLOB EXAMPLES "${PROJECT_DIR}/examples/*.cpp")
foreach(example ${EXAMPLES})
    get_filename_component(name ${example} NAME_WLE)
    add_executable(${name} ${example})
endforeach(example)
