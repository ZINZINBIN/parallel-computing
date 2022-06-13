#include "Control.h"

Control_t Control;
Timer_t Timer;

void Control_t::SetInitialEnv()
{
	omp_set_num_threads(_num_omp_threads);
	MPI_Comm_size(MPI_COMM_WORLD, &_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &_rank);

	stringstream s;
	s << "The number of process  : " << _size;
	Control.Message(s.str()); s.str("");

	s << "The number of threads  : " << _num_omp_threads;
	Control.Message(s.str()); s.str("");

	// Show input file name on terminal
	s << "Input File             : " << Control.GetInputFile();
	Control.Message(s.str()); s.str("");

	s << "Case Id                : " << _caseid;
	Control.Message(s.str()); s.str("");

	_log_fn = _caseid + ".log";
	logFile.open(_log_fn, ios::out);
}

void Control_t::Message(string line, bool is_time, bool is_terminal, bool is_log)
{
	cout.precision(6);
	if (_rank == MASTER_PROCESS) {
		if (is_time) {
			string time_line = Timer.GetCurrentTime();
			if (is_terminal) cout << time_line << "   " << line << endl;
			if (is_log)      logFile << time_line << "   " << line << endl;
		}
		else {
			if (is_terminal) cout << line << endl;
			if (is_log)      logFile << line << endl;
		}
	}
}

void Control_t::Message(string line, int value, bool is_time)
{
	stringstream s;
	s << line << " " << value;
	if (_rank == MASTER_PROCESS) cout << s.str() << endl;
}

void Control_t::Message(string line, double value, bool is_time)
{
	stringstream s;
	s << line << " " << value;
	if (_rank == MASTER_PROCESS) cout << s.str() << endl;
}

void Control_t::Error(ErrorCode error)
{
	switch ((int)error) {
	case (int)ErrorCode::INPUT_OPEN_ERROR:
		Message("Input file does not exists...");
		break;
	case (int)ErrorCode::XSLIB_OPEN_ERROR:
		Message("XS file does not exists...");
		break;
	case (int)ErrorCode::BOOLEAN_ERROR:
		Message("Wrong variable is given in boolean operation...");
		break;
	case (int)ErrorCode::INTEGER_ERROR:
		Message("Wrong variable is given in integer operation...");
		break;
	case (int)ErrorCode::DOUBLE_ERROR:
		Message("Wrong variable is given in double operation...");
		break;
	}

	exit(EXIT_FAILURE);
}

void Control_t::Finalize()
{
	logFile.close();
}

Timer_t::Timer_t()
{
	for (int it = 0; it < (int) Times::NUM_TIMES; it++) {
		_total_times[it] = 0.0;
		_start_times[it] = 0.0;
	}
}

void Timer_t::TimeInit(Times time_type)
{
	_milisecond = 0; _second = 0; _minute = 0; _hour = 0;
	_start_times[(int) time_type] = MPI_Wtime();
}

void Timer_t::TimeStop(Times time_type)
{
	double end_time = MPI_Wtime();
	double elapsed_time = end_time - _start_times[(int) time_type];

	_total_times[(int) time_type]         += elapsed_time;
	_total_times[(int) Times::TOTAL_TIME] += elapsed_time;

	_start_times[(int) time_type] = 0.0;
}

double Timer_t::GetTimeSeconds(Times time_type)
{
	double seconds = _total_times[(int) time_type];

	return seconds;
}

string Timer_t::GetTimeStandard(Times time_type)
{
	int type = (int) time_type;

	_milisecond = _total_times[type] * 1000;

	_second = _milisecond / 1000;
	_milisecond -= 1000 * (_milisecond / 1000);

	_minute = _second / 60;
	_second -= 60 * (_second / 60);

	_hour = _minute / 60;
	_minute -= 60 * (_minute / 60);

	int pad[4];

	pad[0] = _hour < 10 ? 1 : 0;
	pad[1] = _minute < 10 ? 1 : 0;
	pad[2] = _second < 10 ? 1 : 0;
	if (_milisecond < 10)       pad[3] = 2;
	else if (_milisecond < 100) pad[3] = 1;
	else                        pad[3] = 0;

	string time = string(pad[0], '0') + to_string(_hour) + ":"
		+ string(pad[1], '0') + to_string(_minute) + ":"
		+ string(pad[2], '0') + to_string(_second) + "."
		+ string(pad[3], '0') + to_string(_milisecond);

	return time;
}

string Timer_t::GetCurrentTime()
{
	double end_time = MPI_Wtime();
	double total_time = end_time - _start_times[(int) Times::TOTAL_TIME];

	_milisecond = total_time * 1000;

	_second = _milisecond / 1000;
	_milisecond -= 1000 * (_milisecond / 1000);

	_minute = _second / 60;
	_second -= 60 * (_second / 60);

	_hour = _minute / 60;
	_minute -= 60 * (_minute / 60);

	int pad[4];

	pad[0] = _hour < 10 ? 1 : 0;
	pad[1] = _minute < 10 ? 1 : 0;
	pad[2] = _second < 10 ? 1 : 0;
	if (_milisecond < 10)       pad[3] = 2;
	else if (_milisecond < 100) pad[3] = 1;
	else                        pad[3] = 0;

	string time = string(pad[0], '0') + to_string(_hour) + ":"
		+ string(pad[1], '0') + to_string(_minute) + ":"
		+ string(pad[2], '0') + to_string(_second) + "."
		+ string(pad[3], '0') + to_string(_milisecond);

	return time;
}