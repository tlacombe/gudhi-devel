cmake_minimum_required(VERSION 2.6)
project(GUDHI)

find_package(Boost		REQUIRED COMPONENTS)


if(MSVC)
	# Turn off some VC++ warnings
	SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4668 /wd4311 /wd4800 /wd4820 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -std=c++11 -Wall -Wpedantic -lboost_system -lgmpxx -lgmp")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")	
endif()

include_directories(. utils combinatorics geometry topology 
                      Nearest_neighbors/include/gudhi
                      Simplex_tree/include/gudhi
                      Persistent_cohomology/include/gudhi 
                      Hasse_complex/include/gudhi 
                      utilities/include/gudhi )
                      
add_subdirectory(test/skbl/)                   


