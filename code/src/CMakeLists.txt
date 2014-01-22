cmake_minimum_required(VERSION 2.6)
project(GUDHI)

option 		(debug 		"Build GUDHI with debugging on"      	ON)
option 		(optimize 	"Build GUDHI with optimization" 		OFF)

#option        (Boost_NO_SYSTEM_PATHS ON)

#set(Boost_INCLUDE_DIR   ~/ComputationalLibrary/boost_1_53_0/build/include)
#set(Boost_LIBRARY_DIR   ~/$HOME/ComputationalLibrary/boost_1_53_0/build/lib)

find_package(Boost		REQUIRED COMPONENTS)# thread)
#find_package(Doxygen)


add_executable(gudh Filtered_simplex_tree_node.cpp main.cpp)

#Debugging
if 							(debug)
    if 						(optimize)
		set		 			(cxx_flags ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
    else                    (optimize)
            set             (cxx_flags ${CMAKE_CXX_FLAGS_DEBUG})
    endif                   (optimize)
else                        (debug)
    if                      (optimize)
            set             (cxx_flags                  ${CMAKE_CXX_FLAGS_RELEASE})
    else                    (optimize)
            set             (cxx_flags                  ${CMAKE_CXX_FLAGS})
    endif                   (optimize)
endif                       (debug)
add_definitions             (${cxx_flags})

