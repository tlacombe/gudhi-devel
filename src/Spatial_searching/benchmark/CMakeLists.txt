cmake_minimum_required(VERSION 2.6)
project(Spatial_searching_benchmark)

if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()
if (GPROF_PATH)
  # for gprof to make coverage reports - Jenkins
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

# need CGAL 4.8
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.8.0)
    if (EIGEN3_FOUND)      
      add_executable(Spatial_searching_benchmark benchmark_libs.cpp)
      target_link_libraries(Spatial_searching_benchmark 
        ${Boost_DATE_TIME_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
      if (TBB_FOUND)
		add_definitions(-DCGAL_LINKED_WITH_TBB)
        target_link_libraries(Spatial_searching_benchmark ${TBB_LIBRARIES})
      endif(TBB_FOUND)
      
      # NMSLIB
      include_directories("D:/INRIA/Libs/nmslib-1.5.3/similarity_search/include") # CJTODO TEMP
      target_link_libraries(Spatial_searching_benchmark
        optimized "D:/INRIA/Libs/nmslib-1.5.3/similarity_search/x64/Release/MainLib/NonMetricSpaceLib.lib"
        debug "D:/INRIA/Libs/nmslib-1.5.3/similarity_search/x64/Debug/MainLib/NonMetricSpaceLib.lib")
      add_definitions(-DGUDHI_NMSLIB_IS_AVAILABLE)
      
      #FLANN
      # Use the static version ("_s") because the dynamic version is not correctly built on MSVC 
      # (no flann_cpp.lib generated because no exported functions)
      include_directories("D:/INRIA/Libs/flann-4.9.1/src/cpp") # CJTODO TEMP
      target_link_libraries(Spatial_searching_benchmark
        optimized "D:/INRIA/Libs/flann-4.9.1/msvc/lib/Release/flann_cpp_s.lib"
        debug "D:/INRIA/Libs/flann-4.9.1/msvc/lib/Debug/flann_cpp_s.lib")
      add_definitions(-DGUDHI_FLANN_IS_AVAILABLE)
      
      #Nanoflann
      include_directories("D:\INRIA\Libs\nanoflann-1.2.2\include") # CJTODO TEMP
      add_definitions(-DGUDHI_NANOFLANN_IS_AVAILABLE)
      
    endif(EIGEN3_FOUND)
  endif (NOT CGAL_VERSION VERSION_LESS 4.8.0)
endif(CGAL_FOUND)

