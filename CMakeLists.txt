cmake_minimum_required(VERSION 3.10)
project(mhd_solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optional: choose build type if not specified by user
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
endif()

add_subdirectory(math)

add_executable(mhd_solver
    src/euler.cpp
    src/main.cpp
)



# Where to find headers
target_include_directories(mhd_solver PRIVATE
    ${CMAKE_SOURCE_DIR}/include

)




target_compile_options(mhd_solver PRIVATE -Wall -Wextra -Wpedantic -Werror -pg)