cmake_minimum_required(VERSION 3.16)
project(ennesimo)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#boost iostream
find_package(Boost 1.36.0 COMPONENTS system iostreams)
include_directories(${Boost_INCLUDE_DIRS})

#pthread
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_executable(test_gzip test_gzip.cpp)
add_executable(test_unco test_unco.cpp)

target_link_libraries(test_gzip PRIVATE Threads::Threads ${Boost_LIBRARIES})
target_link_libraries(test_unco PRIVATE Threads::Threads)
