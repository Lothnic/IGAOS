enable_testing()

add_executable(linalg_smoke linalg_smoke.cpp)
target_link_libraries(linalg_smoke PRIVATE igaos_simplex)
add_test(NAME linalg_smoke COMMAND linalg_smoke)

add_executable(engine_smoke engine_smoke.cpp)
target_link_libraries(engine_smoke PRIVATE igaos_io igaos_simplex igaos_milp igaos_qp)
target_include_directories(engine_smoke PRIVATE ${CMAKE_SOURCE_DIR}/src/simplex)
add_test(NAME engine_smoke COMMAND engine_smoke)

add_executable(milp_smoke milp_smoke.cpp)
target_link_libraries(milp_smoke PRIVATE igaos_io igaos_simplex igaos_milp)
target_include_directories(milp_smoke PRIVATE ${CMAKE_SOURCE_DIR}/src/simplex)
add_test(NAME milp_smoke COMMAND milp_smoke)

add_executable(warm_smoke warm_smoke.cpp)
target_link_libraries(warm_smoke PRIVATE igaos_io igaos_simplex)
target_include_directories(warm_smoke PRIVATE ${CMAKE_SOURCE_DIR}/src/simplex)
target_compile_definitions(warm_smoke PRIVATE
    IGAOS_SRC_DIR="${CMAKE_SOURCE_DIR}")
add_test(NAME warm_smoke COMMAND warm_smoke)
