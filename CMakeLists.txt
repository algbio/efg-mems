CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(efg-mems CXX)

FIND_PACKAGE(Git QUIET)

cmake_policy(SET CMP0015 NEW)

SET(SDSL_INCLUDE "sdsl-lite/include") #SDSL headeres
SET(SDSL_LIB "sdsl-lite/lib") #SDSL lib

INCLUDE_DIRECTORIES(${SDSL_INCLUDE}) 
LINK_DIRECTORIES(${SDSL_LIB}) 

SET(CMAKE_CXX_STANDARD 11)

SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-g -ggdb -Ofast -fstrict-aliasing -DNDEBUG -march=native")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")

INCLUDE_DIRECTORIES(br-index/src)
INCLUDE_DIRECTORIES(bdbwt/include)

ADD_EXECUTABLE(efg-mems efg-mems.cpp)
TARGET_LINK_LIBRARIES(efg-mems sdsl)
TARGET_LINK_LIBRARIES(efg-mems divsufsort)
TARGET_LINK_LIBRARIES(efg-mems divsufsort64)


