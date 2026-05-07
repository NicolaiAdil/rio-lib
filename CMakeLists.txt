# CMakeLists.txt
cmake_minimum_required(VERSION 3.16)
project(rio_eskf_cpp VERSION 0.1.0 LANGUAGES CXX)

option(RIO_USE_SYSTEM_EIGEN "Prefer system Eigen3 if available" ON)

add_library(rio_eskf_cpp
  src/rio_eskf.cpp
  src/measurements.cpp
)

add_library(rio::rio_eskf_cpp ALIAS rio_eskf_cpp)

target_compile_features(rio_eskf_cpp PUBLIC cxx_std_17)

target_include_directories(rio_eskf_cpp
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Prefer system Eigen on desktop, fallback to vendored headers (important for embedded/submodules).
if(RIO_USE_SYSTEM_EIGEN)
  find_package(Eigen3 QUIET)
endif()

if(Eigen3_FOUND)
  target_include_directories(rio_eskf_cpp PUBLIC ${EIGEN3_INCLUDE_DIRS})
else()
  target_include_directories(rio_eskf_cpp PUBLIC ${CMAKE_CURRENT_LIST_DIR}/third_party/eigen)
endif()

# Consumers may set these via target_compile_definitions() when embedding:
#   RIO_EIGEN_NO_MALLOC
#   RIO_EIGEN_DONT_VECTORIZE
#   RIO_EIGEN_DONT_ALIGN_STATICALLY

include(GNUInstallDirs)

install(TARGETS rio_eskf_cpp
  EXPORT rio_eskf_cppTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT rio_eskf_cppTargets
  FILE rio_eskf_cppTargets.cmake
  NAMESPACE rio::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/rio_eskf_cpp
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${CMAKE_CURRENT_LIST_DIR}/cmake/rio_eskf_cppConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/rio_eskf_cppConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/rio_eskf_cpp
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/rio_eskf_cppConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/rio_eskf_cppConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/rio_eskf_cppConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/rio_eskf_cpp
)