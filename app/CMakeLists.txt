cmake_minimum_required(VERSION 3.10)

project(app VERSION 0.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


####### set optimization levels
if(NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release)
endif()

# set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")



####### main part
add_executable(mysmallpt mysmallpt.cpp)

target_include_directories(mysmallpt
    PRIVATE 
        ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(mysmallpt
    PRIVATE
        smallpt-scene
)
