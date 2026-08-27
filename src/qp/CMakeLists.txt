add_library(igaos_qp STATIC qp.cpp)
target_include_directories(igaos_qp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src/simplex ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(igaos_qp PUBLIC igaos_io)
set_target_properties(igaos_qp PROPERTIES POSITION_INDEPENDENT_CODE ON)
