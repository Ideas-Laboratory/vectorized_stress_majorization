
set(${CMAKE_CXX_FLAGS} "-D_FORCE_INLINES ${CMAKE_CXX_FLAGS}")
project(test)

cmake_minimum_required(VERSION 2.8.5)

find_package(CUDA)
#enable_language(CUDA)
find_package(GLUT)

include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )


#if(NOT CMAKE_BUILD_TYPE)
#	set(CMAKE_BUILD_TYPE Release Debug FORCE)
#	message(STATUS "build type not set ${CMAKE_BUILD_TYPE}")
#endif()
#add_executable(test1 test.c)
#set_source_files_properties(test.c PROPERTIES LANGUAGE CXX )
list(APPEND CUDA_NVCC_FLAGS "-arch=sm_35;-std=c++11;-O3;-DVERBOSE;-D_FORCE_INLINES;-D_MWAITXINTRIN_H_INCLUDED")
#message(STATUS "cuda flags: ${CMAKE_CUDA_FLAGS}")
#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_35 -std=c++11 -DVERBOSE -D_FORCE_INLINES -D_MWAITXINTRIN_H_INCLUDED")
#message(STATUS "cuda flags modified: ${CMAKE_CUDA_FLAGS}")
#message(STATUS "${cuda_target}")

#add_library(${cuda_target})

CUDA_INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/")

CUDA_INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/Eigen")
SET(CUDA_PROPAGATE_HOST_FLAGS OFF)
SET(CUDA_SEPARABLE_COMPILATION ON)
SET(CUDA_HOST_COMPILER g++)
#SET(CUDA_HOST_COMPILER cl)
file(GLOB CUDA_FILES "${CMAKE_CURRENT_SOURCE_DIR}/" *.cu *.c *.cpp)
#file(GLOB CXX_FILES "${CMAKE_CURRENT_SOURCE_DIR}/"

CUDA_COMPILE(CU_O ${CUDA_FILES})
CUDA_ADD_LIBRARY(KERNELS ${CUDA_FILES} STATIC)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/Eigen/")
message(STATUS "${CMAKE_CURRENT_SOURCE_DIR}/Eigen/")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/")
include_directories("${CMAKE_CURRENT_DIR}/")
include_directories("${CUDA_SDK_ROOT_DIR}/include/")
#LINK_DIRECTORIES("${CUDA_SDK_ROOT_DIR}/lib/")
#add_library(cudart.a)
cuda_add_executable(test1 ${CUDA_FILES} OPTIONS -gencode arch=compute_35,code=sm_35)
#add_executable(test1 ${CUDA_FILES})
target_link_libraries(test1 ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )
#target_link_libraries(test1 -lGL)
target_link_libraries(test1 ${CUDA_LIBRARIES} -lGL)

message(STATUS "${GLUT_LIBRARY}")
message(STATUS "cuda files ${CUDA_FILES}")
message(STATUS "${CUDA_NVCC_FLAGS}")
#set(CMAKE_EXE_LINKER_FLAGS -lcudart -lcudadevrt)
