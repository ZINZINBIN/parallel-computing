#pragma once
// Headers
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <cmath>
#ifdef __linux__
#include <string.h>
#endif
/// Parallel
#include <omp.h>
#include <mpi.h>

using namespace std;


// Define class and struct
/// Util.h
template <typename T> class Array_t;
template <typename T> class CSR_t;
/// Geometry.h
struct Assembly_t;

// Declare base variables
const char  BLANK  = ' ';
const char  PANG   = '!';  /// For annotation

// Define Parameters
#define INVALID        -1
#define MASTER_PROCESS 0 /// Define master process of mpi
/// Direction
#define WEST  0
#define EAST  1
#define NORTH 2
#define SOUTH 3
#define LOWER 4
#define UPPER 5

// Define function
#ifdef _WIN32            /// Enviroment set up 
#define strtok strtok_s
#else
#define strtok strtok_r
#endif

// Define enum class
/// Util.h
enum class ErrorCode {
	INPUT_OPEN_ERROR,
	XSLIB_OPEN_ERROR,
	BOOLEAN_ERROR,
	INTEGER_ERROR,
	DOUBLE_ERROR,
};
/// Control.h 
enum class SolverType {
	BICGSTAB,
	SOR,
};
enum class Times {
	GEOMETRY_TIME,
	SOLVER_SETTING_TIME,
	SOLVER_TIME,
	TOTAL_TIME,
	NUM_TIMES,
};