enable_testing()

add_executable(linalg_smoke linalg_smoke.cpp)
target_link_libraries(linalg_smoke PRIVATE igaos_simplex)
add_test(NAME linalg_smoke COMMAND linalg_smoke)

add_executable(engine_smoke engine_smoke.cpp)
target_link_libraries(engine_smoke PRIVATE igaos_io igaos_simplex igaos_milp igaos_qp)
target_include_directories(engine_smoke PRIVATE ${CMAKE_SOURCE_DIR}/src/simplex)
add_test(NAME engine_smoke COMMAND engine_smoke)
