cmake_minimum_required(VERSION 3.24)
project(igaos VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Default to Release: an unoptimized build (empty CMAKE_BUILD_TYPE) left the
# dense LU/eta kernels at -O0 — ~9x slower, and the Netlib 30s gate results
# depended on it silently. No build type was ever set anywhere, so every
# default configure was running unoptimized.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include(CheckLanguage)
# IGAOS_ENABLE_CUDA: AUTO detects the toolkit (default), ON requires it,
# OFF builds CPU-only — used by pip config-settings for source builds.
set(IGAOS_ENABLE_CUDA "AUTO" CACHE STRING "CUDA engine: AUTO|ON|OFF")
set(IGAOS_HAVE_PDHG 0)
if(NOT IGAOS_ENABLE_CUDA STREQUAL "OFF")
  check_language(CUDA)
  if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    set(IGAOS_HAVE_PDHG 1)
  elseif(IGAOS_ENABLE_CUDA STREQUAL "ON")
    message(FATAL_ERROR "IGAOS_ENABLE_CUDA=ON but no CUDA toolchain found")
  else()
    message(STATUS "CUDA toolchain not found: building CPU-only")
  endif()
else()
  message(STATUS "IGAOS_ENABLE_CUDA=OFF: building CPU-only")
endif()

add_subdirectory(src/common)
add_subdirectory(src/linalg)
add_subdirectory(src/simplex)
add_subdirectory(src/milp)
add_subdirectory(src/qp)
add_subdirectory(src/io)
if(IGAOS_HAVE_PDHG)
  add_subdirectory(src/pdhg)
endif()
add_subdirectory(src/api)

# wheel builds (scikit-build-core) need only the python module — the CLI
# and test executables are source-build artifacts
if(NOT SKBUILD)
    enable_testing()
    add_subdirectory(tests)
endif()
