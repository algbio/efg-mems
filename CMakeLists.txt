CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(efg-mems CXX)

FIND_PACKAGE(Git QUIET)

SET(SDSL_INCLUDE "~/include") #SDSL headeres
SET(SDSL_LIB "~/lib") #SDSL lib

INCLUDE_DIRECTORIES(${SDSL_INCLUDE}) 
LINK_DIRECTORIES(${SDSL_LIB}) 

SET(CMAKE_CXX_STANDARD 11)

SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-g -ggdb -Ofast -fstrict-aliasing -DNDEBUG -march=native")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")

INCLUDE_DIRECTORIES(../br-index/src)

ADD_EXECUTABLE(efg-mems efg-mems.cpp)
TARGET_LINK_LIBRARIES(efg-mems sdsl)
TARGET_LINK_LIBRARIES(efg-mems divsufsort)
TARGET_LINK_LIBRARIES(efg-mems divsufsort64)


