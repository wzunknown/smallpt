cmake_minimum_required(VERSION 3.10)

project(mysmallpt VERSION 0.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


####### set optimization levels
if(NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


####### find dependencies
find_package(OpenMP REQUIRED)
# find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# find_package(yaml-cpp)

add_subdirectory(lib)
# add_subdirectory(src)



####### main part
file(GLOB SRC_FILE src/*.cpp)
# set(SRC_FILE src/color.cpp src/scene.cpp src/sphere.cpp src/vector.cpp)

add_executable(${PROJECT_NAME} app/mysmallpt.cpp ${SRC_FILE})

target_include_directories(${PROJECT_NAME} 
    PRIVATE 
        include
)

target_link_libraries(${PROJECT_NAME}
    PRIVATE
        yaml-cpp
)