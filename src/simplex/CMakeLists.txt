add_library(igaos_simplex STATIC simplex.cpp)
target_include_directories(igaos_simplex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(igaos_simplex PUBLIC igaos_io)
