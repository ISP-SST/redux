# Global configuration for Apple/Mac
#
# Author: Tomas Hillberg

MESSAGE(STATUS "Loading mac-specific configuration. (${CMAKE_CURRENT_LIST_FILE})")

# TODO would be neater if system defaults worked without this hardcoding
set(CMAKE_FIND_LIBRARY_SUFFIXES ".dylib" ".a")
