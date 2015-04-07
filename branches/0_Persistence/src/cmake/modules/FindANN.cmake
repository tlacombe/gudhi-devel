# - Try to find the ANN libraries
# This module defines:
#  ANN_FOUND             - system has ANN lib
#  ANN_INCLUDE_DIR       - the ANN include directory
#  ANN_LIBRARIES_DIR     - directory where the ANN libraries are located
#  ANN_LIBRARIES         - Link these to use ANN

# TODO: support MacOSX

include(FindPackageHandleStandardArgs)

if(ANN_INCLUDE_DIR)
  set(ANN_in_cache TRUE)
else()
  set(ANN_in_cache FALSE)
endif()
if(NOT ANN_LIBRARIES)
  set(ANN_in_cache FALSE)
endif()

# Is it already configured?
if (ANN_in_cache)
  set(ANN_FOUND TRUE)
else()
  find_path(ANN_INCLUDE_DIR
            NAMES ANN/ANN.h
            HINTS ENV ANN_INC_DIR
                  ENV ANN_DIR
            PATH_SUFFIXES include
  	        DOC "The directory containing the ANN header files"
           )

  find_library(ANN_LIBRARIES NAMES ANN libANN
    HINTS ENV ANN_LIB_DIR
          ENV ANN_DIR
    PATH_SUFFIXES lib
    DOC "Path to the ANN library"
    )

  if ( ANN_LIBRARIES )
    get_filename_component(ANN_LIBRARIES_DIR ${ANN_LIBRARIES} PATH CACHE )
  endif()

  # Attempt to load a user-defined configuration for ANN if couldn't be found
  if ( NOT ANN_INCLUDE_DIR OR NOT ANN_LIBRARIES_DIR )
    include( ANNConfig OPTIONAL )
  endif()

  find_package_handle_standard_args(ANN "DEFAULT_MSG" ANN_LIBRARIES ANN_INCLUDE_DIR)

endif()
