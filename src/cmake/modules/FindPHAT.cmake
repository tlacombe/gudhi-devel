# - Try to find the PHAT libraries
# This module defines:
#  PHAT_FOUND             - system has PHAT lib
#  PHAT_INCLUDE_DIR       - the PHAT include directory


include(FindPackageHandleStandardArgs)

if(PHAT_INCLUDE_DIR)
  set(PHAT_in_cache TRUE)
else()
  set(PHAT_in_cache FALSE)
endif()

# Is it already configured?
if (PHAT_in_cache)

  set(PHAT_FOUND TRUE)

else()

  find_path(PHAT_INCLUDE_DIR
            NAMES phat/boundary_matrix.h
            HINTS ENV PHAT_INC_DIR
                  ENV PHAT_DIR
            PATH_SUFFIXES include
  	        DOC "The directory containing the PHAT header files"
           )

  find_package_handle_standard_args(PHAT "DEFAULT_MSG" PHAT_INCLUDE_DIR)

endif()
