add_library(igaos_milp STATIC milp.cpp)
target_include_directories(igaos_milp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(igaos_milp PUBLIC igaos_simplex)
set_target_properties(igaos_milp PROPERTIES LINKER_LANGUAGE CXX)
