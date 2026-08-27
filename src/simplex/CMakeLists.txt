add_library(igaos_simplex STATIC simplex.cpp)
set_target_properties(igaos_simplex PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(igaos_simplex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(igaos_simplex PUBLIC igaos_io)
