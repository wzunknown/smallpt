cmake_minimum_required(VERSION 3.10)

project(mysmallpt-mpi VERSION 0.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


####### set optimization levels
if(NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release)
endif()

# set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

####### check dependencies
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

####### main part
add_executable(mysmallpt-mpi mysmallpt-mpi.cpp)


target_include_directories(mysmallpt-mpi
    PRIVATE 
        ${CMAKE_SOURCE_DIR}/include
        MPI_CXX
)

target_link_libraries(mysmallpt-mpi
    PRIVATE
        smallpt-scene
)
