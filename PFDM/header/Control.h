#pragma once
#include "Defines.h"

/// Control.h
class Control_t;
class Timer_t;

extern Control_t Control;
extern Timer_t Timer;

class Control_t
{
private:
	int _rank = 0, _size = 1; /// MPI Parameters
	int _num_omp_threads = 1;
	int _max_outer, _max_inner;
	double _eigen_crit, _psi_crit, _inner_crit;
	string _caseid;
	string _input_fn, _log_fn, _xslib_fn;
	ofstream logFile;
	SolverType _solver_type = SolverType::BICGSTAB;

public:
	void ReadCMDLine(int argc, char** argv);
	void ScanNumThread(ifstream& file);
	void SetInitialEnv();
	void Message(string line, bool is_time = true, bool is_terminal = true, bool is_log = true);
	void Message(string line, int value, bool is_time = false);
	void Message(string line, double value, bool is_time = false);
	void Error(ErrorCode error);
	void Finalize();
	// Setting functions
	void SetInputFile(string fn) { _input_fn = fn; };
	void SetXSFile(string fn)    { _xslib_fn = fn; };
	void SetCaseID(string field) { _caseid = field; };
	/// Set convergence criteria
	void SetMaxOuter(int max_outer)            { _max_outer = max_outer; };
	void SetMaxInner(int max_inner)            { _max_inner = max_inner; };
	void SetEigenCriteria(double eigen_crit)   { _eigen_crit = eigen_crit; };
	void SetPsiCriteria(double psi_crit)       { _psi_crit = psi_crit; };
	void SetInnerCriteria(double inner_crit)   { _inner_crit = inner_crit; };
	void SetSolverType(SolverType solver_type) { _solver_type = solver_type; };
	// Query functions
	int GetMPIRank()             { return _rank; };
	int GetMPISize()             { return _size; };
	double GetPsiCritera()       { return _psi_crit; };
	double GetEigenCritera()     { return _eigen_crit; };
	double GetInnerCriteria()    { return _inner_crit; };
	string GetCaseID()           { return _caseid; };
	string GetInputFile()        { return _input_fn; };
	string GetXSFile()           { return _xslib_fn; };
	SolverType GetSolverType()   { return _solver_type;};
};

class Timer_t
{
private:
	int _milisecond, _second, _minute, _hour;
	double _total_times[(int) Times::NUM_TIMES];
	double _start_times[(int) Times::NUM_TIMES];

public:
	Timer_t();
	void Initialize() {_start_times[(int) Times::TOTAL_TIME] = MPI_Wtime(); };
	void TimeInit(Times time_type);
	void TimeStop(Times time_type);
	double GetTimeSeconds(Times time_type);
	string GetTimeStandard(Times time_type);
	string GetCurrentTime();
};