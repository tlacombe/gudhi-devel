cmake_minimum_required(VERSION 2.6)
project(GUDHIskbl)

find_package(Boost)

if(MSVC)
	# Turn off some VC++ warnings
	SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4668 /wd4311 /wd4800 /wd4820 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -std=c++11 -Wall -lboost_system -lgmpxx -lgmp")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

if(NOT Boost_FOUND)
  message(FATAL_ERROR "NOTICE: This demo requires Boost and will not be compiled.")
else()
	INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
	LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
endif()

include_directories(../../ ../../utils ../../combinatorics ../../geometry ../../topology)

add_executable ( TestSkeletonBlockerComplex TestSkeletonBlockerComplex.cxx )
add_executable ( TestSimplifiable           TestSimplifiable.cxx )
add_executable ( TestContraction            TestContraction.cxx )
