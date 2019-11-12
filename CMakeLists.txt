cmake_minimum_required(VERSION 2.6)
project(thetastarcpoc)

include(ExternalProject)
find_package(PkgConfig)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/Modules)

find_package(PROJ4 REQUIRED)
find_package(GEOS REQUIRED)

if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
    ADD_COMPILE_OPTIONS(-std=c++17 -g ${GST_CFLAGS})
endif()

SET(THREADS_PREFER_PTHREAD_FLAG ON)


add_executable(thetastarcpoc main.cpp)

include_directories($include_directories ${GEOS_INCLUDE_DIR} ${PROJ4_INCLUDE_DIRS})

TARGET_LINK_LIBRARIES(thetastarcpoc ${GEOS_LIBRARY} ${PROJ4_LIBRARY})


install(TARGETS thetastarcpoc RUNTIME DESTINATION bin)


