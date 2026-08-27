add_executable(igaos igaos_cli.cpp)
target_link_libraries(igaos PRIVATE igaos_io igaos_simplex igaos_milp)
target_include_directories(igaos PRIVATE ${CMAKE_SOURCE_DIR}/include)

if(IGAOS_HAVE_PDHG)
    target_link_libraries(igaos PRIVATE igaos_pdhg)
    target_compile_definitions(igaos PRIVATE IGAOS_HAS_PDHG)
endif()

# Python bindings (opt-in; dev tool pybind11 is allowed per docs/DEPENDENCIES.md)
find_package(pybind11 CONFIG QUIET)
if(pybind11_FOUND)
    # Target name differs from the CLI executable; OUTPUT_NAME makes the
    # built file import as `igaos` (matches PYBIND11_MODULE name).
    pybind11_add_module(igaos_py ${CMAKE_SOURCE_DIR}/python/igaos_bindings.cpp)
    target_link_libraries(igaos_py PRIVATE igaos_io igaos_simplex igaos_milp)
    target_include_directories(igaos_py PRIVATE ${CMAKE_SOURCE_DIR}/src/api)
    if(IGAOS_HAVE_PDHG)
        target_link_libraries(igaos_py PRIVATE igaos_pdhg)
        target_compile_definitions(igaos_py PRIVATE IGAOS_HAS_PDHG)
    endif()
    set_target_properties(igaos_py PROPERTIES
        OUTPUT_NAME igaos
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python)
    message(STATUS "pybind11 found: building python/igaos bindings")
else()
    message(STATUS "pybind11 not found: skipping python module")
endif()
