cmake_minimum_required(VERSION 3.10)

project(smallpt-scene VERSION 0.1)


option(USE_LOCAL_YAML_CPP "use user-installed yaml-cpp library" OFF)

###### --std=c++11
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
file(GLOB SMALLPT_SRC ${PROJECT_SOURCE_DIR}/src/*.cpp)
# set(SMALLPT_SRC src/color.cpp src/scene.cpp src/sphere.cpp src/vector.cpp)

configure_file(${CMAKE_SOURCE_DIR}/default.yaml ${CMAKE_BINARY_DIR}/default.yaml COPYONLY)

add_library(smallpt-scene STATIC ${SMALLPT_SRC})

target_include_directories(smallpt-scene
    PRIVATE 
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(smallpt-scene
    PRIVATE
        yaml-cpp
)

####### find dependencies
# find_package(OpenMP)
find_package(OpenMP COMPONENTS CXX)
if (OPENMP_FOUND)
    target_link_libraries(smallpt-scene
        PRIVATE
            OpenMP::OpenMP_CXX
    )
else()
    message(STATUS "Not found OpenMP for cxx, program will be slow")
endif()
# if (OPENMP_FOUND)
#     set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()

if (USE_LOCAL_YAML_CPP)
    find_library(YAML_CPP_LIBRARY NAMES yaml-cpp)
endif()
if (YAML_CPP_LIBRARY)
    message(STATUS "Found libyaml-cpp")
else()
    message(STATUS "Not found yaml-cpp, will build from source in lib/yaml-cpp")
    add_subdirectory(${PROJECT_SOURCE_DIR}/lib)
endif()

add_subdirectory(app)


