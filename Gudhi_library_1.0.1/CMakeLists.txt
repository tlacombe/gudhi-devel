cmake_minimum_required(VERSION 2.6)
project(Gudhi)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules)

SET(BOOST_MIN_VERSION "1.48.0")
find_package(Boost ${BOOST_MIN_VERSION}	REQUIRED COMPONENTS program_options system )
find_package(GMPXX        QUIET)

if (NOT CMAKE_BUILD_TYPE)
 set(CMAKE_BUILD_TYPE "Release")
endif()

include_directories( .
		     ${Boost_INCLUDE_DIRS}
		     ${GMPXX_INCLUDE_DIR}
                     Simplex_tree/include/gudhi
                     Persistent_cohomology/include/gudhi 
                     Hasse_complex/include/gudhi 
                     utilities/include/gudhi )

add_executable( rips_persistence examples/rips_persistence.cpp )
target_link_libraries( rips_persistence ${Boost_LIBRARIES} )

if(MSVC)
else()
  list(APPEND CMAKE_CXX_FLAGS "-std=c++0x")
endif()

if( GMP_FOUND ) 

    add_executable ( rips_multifield_persistence examples/rips_multifield_persistence.cpp )
    target_link_libraries( rips_multifield_persistence ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} )

    add_executable( performance_rips_persistence examples/performance_rips_persistence.cpp )
    target_link_libraries( performance_rips_persistence ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} )

else()

endif()
