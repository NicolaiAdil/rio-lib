#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "rio::rio_eskf_cpp" for configuration ""
set_property(TARGET rio::rio_eskf_cpp APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(rio::rio_eskf_cpp PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/librio_eskf_cpp.a"
  )

list(APPEND _cmake_import_check_targets rio::rio_eskf_cpp )
list(APPEND _cmake_import_check_files_for_rio::rio_eskf_cpp "${_IMPORT_PREFIX}/lib/librio_eskf_cpp.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
