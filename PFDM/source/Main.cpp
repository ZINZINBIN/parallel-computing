#include "Defines.h"
#include "Control.h"
#include "IO.h"
#include "Geometry.h"
#include "Solver.h"

void Init()
{
	Timer.Initialize();

	Control.SetInitialEnv();
}

void Solver()
{
	if (Control.GetSolverType() == SolverType::BICGSTAB) {

		BiCGSTAB.Setting();
		BiCGSTAB.Solve();
		 
	}
	else if (Control.GetSolverType() == SolverType::SOR) {
		SOR.Setting();
		SOR.Solve();
	}
}

void Finalize()
{
	PrintOutput();

	MPI_Finalize();
	Control.Message("Calculation is finished...", false);

	Control.Message("====================================================================================", false);
	stringstream s;
	if (Control.GetSolverType() == SolverType::BICGSTAB) 
		s << "                                     keff  =  " <<  BiCGSTAB.GetKeffective();
	else  
		s << "                                     keff  =  " << BiCGSTAB.GetKeffective();
	Control.Message(s.str(), false); s.str("");
	Control.Message("====================================================================================", false);
	s << "                             - Total Time : " << Timer.GetTimeSeconds(Times::TOTAL_TIME); Control.Message(s.str(), false); s.str("");
	s << "                               - Geometry Time       : " << Timer.GetTimeSeconds(Times::GEOMETRY_TIME); Control.Message(s.str(), false); s.str("");
	s << "                               - Solver Setting Time : " << Timer.GetTimeSeconds(Times::SOLVER_SETTING_TIME); Control.Message(s.str(), false); s.str("");
	s << "                               - Solver Time         : " << Timer.GetTimeSeconds(Times::SOLVER_TIME); Control.Message(s.str(), false); s.str("");

	Control.Finalize();
}

int main(int argv, char* argc[])
{
	// Initialize enviornment
	Control.ReadCMDLine(argv, argc);

	MPI_Init(&argv, &argc);

	Init();

	// Read input file
	ReadInput();

	// Read XS library
	ReadXSLib();

	// Geometry configuration
	Geometry.ConstructGeometry();
	
	// Solver
	Solver();
	
	//// Finalize program
	Finalize();

	return EXIT_SUCCESS;
}