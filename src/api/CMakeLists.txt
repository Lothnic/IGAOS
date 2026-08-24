add_executable(igaos igaos_cli.cpp)
target_link_libraries(igaos PRIVATE igaos_io)
target_include_directories(igaos PRIVATE ${CMAKE_SOURCE_DIR}/include)
