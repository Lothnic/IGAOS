enable_testing()

add_executable(linalg_smoke linalg_smoke.cpp)
target_link_libraries(linalg_smoke PRIVATE igaos_simplex)
add_test(NAME linalg_smoke COMMAND linalg_smoke)
