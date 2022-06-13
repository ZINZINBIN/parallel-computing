#pragma once
#include "Defines.h"
#include "Util.h"
#include "Control.h"
#include "Geometry.h"
#include "Xsec.h"

/// Solver.h
class Solver_t;
class BiCGSTAB_t;
class SOR_t;

extern class BiCGSTAB_t BiCGSTAB;
extern class SOR_t SOR;

class Solver_t
{
protected:
	// Geometry Informations
	double _volume;
	double _surface_size[3];
	double _boundary_diffusivity[6];
	Array_t<int> _mesh_materials;
	Array_t<int> _adjacent_nodes;
	Array_t<double> _diffusivity;
	
	// Solver Informations
	int _unit_num_nodes;
	int _num_groups, _num_nodes, _node_offset;
	int _num_rows, _start_row, _end_row, _row_offset;
	int _unit_buffer_size;
	int _buffer_nodes_size;
	// Parameters
	bool _is_convergence = false;
	double _keff, _prev_keff;
	CSR_t<double> _M;
	Array_t<double> _chi_arr;
	Array_t<double> _psi, _prev_psi;
	Array_t<double> _phi, _prev_phi;
	Array_t<double> _RHS;
	Array_t<double> _total_phi, _assembly_phi;

	void GeometryInfoSetting();
	void MPISetting();
	void InitializeBasicArray();
	void ConstructMatrix();
	void UpdateRHS(int icycle);
	void UpdateFissionSource(int icycle);
	void DetermineConvergence(int icycle);
	void GatherTotalFlux();
	void GatherAssemblyFlux();
	virtual void Initialize()   = 0;
	virtual void LinearSolver(int icycle) = 0;

public:
	void Setting();
	void Solve();
	// Query Functions
	double GetKeffective() { return _keff; };
	double GetAssemblyFlux(int ig, int iasyX, int iasyY, int iasyZ) { return _assembly_phi(ig, iasyX, iasyY, iasyZ); };
};

class BiCGSTAB_t: public Solver_t
{
private:
	double _output_level = 1;
	double _zero_residual_size;
	double _rho, _prev_rho;
	double _omega, _prev_omega;
	double _alpha;
	double _beta;
	double _temp_value_buffer;
	bool _is_inner_convergence     = false;
	Array_t<double> _temp_arr_buffer;
	Array_t<double> _zero_residual, _residual, _prev_residual, _hat_residual;
	Array_t<double> _nu, _prev_nu;
	Array_t<double> _pvalue, _prev_pvalue;
	Array_t<double> _hash, _svalue, _tvalue;

	virtual void Initialize();
	void AllocateArray();
	void UpdatePvalue();
	void UpdateHash();
	void UpdateResidual(int icycle);

public:
	void LinearSolver(int icycle);
};

/**
 * @brief Successive Over-relaxation method
 * @date 2022-05-30
 */

class SOR_t: public Solver_t
{
	private:
		double _w_optimal;
		double _spectral_radius;
		double _beta;
		double _temp_value_buffer;
		double _diff_norm;
		double _eps;

		bool _is_inner_convergence = false;
		int _max_inner_iterations; 

		Array_t<double> _matmul_buffer; // sum of a_ij * phi_j 
		Array_t<double> _residual_buffer; // b - sum of a_ij * phi_j

		virtual void Initialize(); // initialize psi, psi_prev, matrix A, vector B, etc..
		void AllocateArray(); // allocate array
		void CalculateWvalue(); // calculate optimal w from matrix A only once
		void UpdateVector(); // update phi using SOR algorithm per 1 epoch
		void CheckIsConvergence(); // Check if _phi_prev and _phi are converged
		double CalculateDiffNorm(); // calculate L2 Norm between psi_prev and psi for convergence

	public:
		void LinearSolver(int icycle); // update psi from psi_prev for 1 epoch
};