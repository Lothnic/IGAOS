add_library(igaos_pdhg STATIC pdhg.cu)
target_include_directories(igaos_pdhg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(igaos_pdhg PUBLIC igaos_io CUDA::cudart CUDA::cublas CUDA::cusparse)
set_target_properties(igaos_pdhg PROPERTIES CUDA_ARCHITECTURES native)
