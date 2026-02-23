# rio-lib (rio_eskf_cpp)

Fixed-size Eigen-based ESKF for radar-inertial odometry.

## Usage

Add as a git submodule and in CMake:

```cmake
add_subdirectory(rio-lib)
target_link_libraries(your_target PRIVATE rio::rio_eskf_cpp)