cmake_minimum_required(VERSION 2.6)
project(Simplex_array_list_tests)

add_executable ( SalUT test.cpp )
target_link_libraries(SalUT ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})


# Unitary tests
add_test(NAME SalUT 
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/SalUT
    ${CMAKE_SOURCE_DIR}/src/Simplex_array_list/test/test.txt
    # XML format for Jenkins xUnit plugin 
    --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/SalUT.xml --log_level=test_suite --report_level=no)

