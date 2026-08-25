add_executable(igaos igaos_cli.cpp)
target_link_libraries(igaos PRIVATE igaos_io igaos_simplex igaos_milp)
target_include_directories(igaos PRIVATE ${CMAKE_SOURCE_DIR}/include)

if(IGAOS_HAVE_PDHG)
    target_link_libraries(igaos PRIVATE igaos_pdhg)
    target_compile_definitions(igaos PRIVATE IGAOS_HAS_PDHG)
endif()
