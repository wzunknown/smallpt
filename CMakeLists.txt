cmake_minimum_required(VERSION 3.10)

project(mysmallpt VERSION 0.1)

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

####### find dependencies
# find_package(OpenMP)
find_package(OpenMP COMPONENTS CXX)
if (OPENMP_FOUND)
    target_link_libraries(${PROJECT_NAME}
    PRIVATE
    OpenMP::OpenMP_CXX
    )
else()
    message(STATUS "Not found OpenMP for cxx, process will be slow")
endif()
# if (OPENMP_FOUND)
#     set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()

find_library(YAML_CPP_LIBRARY NAMES yaml-cpp)
# find_library(LIB)
if (YAML_CPP_LIBRARY)
    message(STATUS "Found libyaml-cpp")
else()
    message(STATUS "Not found yaml-cpp, will build from source in lib/yaml-cpp")
    add_subdirectory(lib)
endif()

# add_subdirectory(src)


