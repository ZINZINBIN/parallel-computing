### PFDM CMake Script

set (CMAKE_C_COMPILER mpicc)
set (CMAKE_CXX_COMPILER mpic++)
set (CMAKE_CXX_STANDARD 17)

project (
	PFDM
)

find_package(OpenMP REQUIRED)

set (CC -gencode arch=compute_75,code=sm_75)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O2 -w  ${OpenMP_CXX_FLAGS}")

include_directories (
	./source
	./header
)

file (GLOB CPP_SOURCE_FILES ./source/*.cpp)
file (GLOB CPP_HEADER_FILES ./header/*.h)

add_executable (
	PFDM
	${CPP_SOURCE_FILES}
	${CPP_HEADER_FILES}
)

target_link_libraries (
	PFDM
	stdc++fs
)
