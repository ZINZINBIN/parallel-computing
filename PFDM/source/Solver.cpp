#include "Solver.h"

BiCGSTAB_t BiCGSTAB;
SOR_t SOR;

void Solver_t::Setting()
{
	Timer.TimeInit(Times::SOLVER_SETTING_TIME);

	Control.Message("BiCGSTAB Method is adopted...");

	GeometryInfoSetting();

	MPISetting();

	InitializeBasicArray();

	ConstructMatrix();

	Timer.TimeStop(Times::SOLVER_SETTING_TIME);
}

void Solver_t::GeometryInfoSetting()
{
	_volume = Geometry.GetUnitVolume();
	Geometry.GetMeshMaterial(_mesh_materials);
	Geometry.GetDiffusivity(_diffusivity);
	Geometry.GetSurfaceSize(_surface_size);
	Geometry.GetAdjacentNode(_adjacent_nodes);
	Geometry.GetBoundaryDiffusivity(_boundary_diffusivity);
}

void Solver_t::MPISetting()
{
	_unit_num_nodes      = Geometry.GetNumMeshPerFloor();
	_num_groups          = XS.GetNumOfGruop();

	_num_nodes           = _unit_num_nodes * Geometry.GetNumLayer();
	_node_offset         = Geometry.GetStartLayer() * _unit_num_nodes;
	
	_unit_buffer_size    = _unit_num_nodes * _num_groups;
	_buffer_nodes_size   = _num_nodes + 2 * _unit_num_nodes;

	_num_rows            = _unit_buffer_size * Geometry.GetNumLayer();
	_start_row           = Geometry.GetStartLayer() * _unit_buffer_size;
	_end_row             = _start_row + _num_rows - 1;
	
	if (Control.GetMPIRank() == MASTER_PROCESS) _row_offset = 0;
	else                                        _row_offset = _unit_buffer_size;

	if (Control.GetMPIRank() == MASTER_PROCESS)                _buffer_nodes_size -= _unit_num_nodes;
	else if (Control.GetMPIRank() == Control.GetMPISize() - 1) _buffer_nodes_size -= _unit_num_nodes;
}


void Solver_t::InitializeBasicArray()
{
	Control.Message("Initialize array...");
	_keff = 1.0; _prev_keff = 1.0;
	_phi.Create(_num_groups, _buffer_nodes_size);       _phi = 0.0; 
	for (int irow = 0; irow < _num_rows; irow++) _phi(irow + _row_offset) = 1.0;
	_prev_phi.Create(_num_groups, _buffer_nodes_size);  _prev_phi = 1.0;
	
	_psi.Create(_num_nodes); _prev_psi.Create(_num_nodes);

	_chi_arr.Create(_num_groups, _num_nodes);           _chi_arr = 0.0;

	_M.SetRowSize(_num_rows);
	_RHS.Create(_num_groups, _num_nodes);               _RHS = 0.0;

	double volume = Geometry.GetUnitVolume();
	#pragma omp parallel for schedule(guided)
	for (int inode = 0; inode < _num_nodes; inode++) {
		int inode_global = inode + _node_offset;
		int imat = _mesh_materials(inode_global);
		Mat_t material = XS.GetMaterialInfo(imat);

		double sum = 0.0;
		for (int ig = 0; ig < _num_groups; ig++) {
			int idx = inode * _num_groups + ig + _row_offset;
			sum += material._nusigf[ig] * _phi(idx);
		}

		_psi(inode) = sum * volume;
	}

	_prev_psi = _psi;
}

void Solver_t::ConstructMatrix()
{
	Control.Message("Construct matrix...");
	int adjacent_node_idx[6];
	double volume = Geometry.GetUnitVolume();
	double dTilda[6];

	for (int inode = 0; inode < _num_nodes; inode++) {
		int inode_global = inode + _node_offset;
		int imat         = _mesh_materials(inode_global);
		Mat_t material   = XS.GetMaterialInfo(imat);

		for (int ig = 0; ig < _num_groups; ig++) {
			double diagonal_element = 0.0;
			int matrix_idx          = inode_global * _num_groups + ig;
			int matrix_idx_local    = inode * _num_groups + ig;

			for (int idim = 0; idim < 3; idim++) {
				double surface = _surface_size[idim];
				double betaCur = _diffusivity(idim, ig, inode_global);
				double betaAdj;

				for (int idir = 0; idir < 2; idir++) {
					int iadj    = 2 * idim + idir;
					int adj_idx = _adjacent_nodes(iadj, inode_global);

					if (adj_idx == INVALID) {
						adjacent_node_idx[iadj] = adj_idx;
						betaAdj = _boundary_diffusivity[iadj];
					}
					else {
						adjacent_node_idx[iadj] = adj_idx * _num_groups + ig;
						betaAdj = _diffusivity(idim, ig, adj_idx);
					}

					dTilda[iadj]      = 2 * surface * betaAdj * betaCur / (betaAdj + betaCur);
					diagonal_element += dTilda[iadj];
				}
			}

			diagonal_element += material._sigr[ig] * volume;

			for (int idx = 4; idx >= 0; idx -= 2) {
				if (adjacent_node_idx[idx] == INVALID) continue;
				_M.PushBack(matrix_idx_local, adjacent_node_idx[idx], -dTilda[idx]);
			}

			for (int igg = 0; igg < _num_groups; igg++) {
				int group_matrix_idx = inode_global * _num_groups + igg;
				double temp_scat     = material._scat(igg, ig);

				if (igg == ig) {
					_M.PushBack(matrix_idx_local, group_matrix_idx, diagonal_element);
				}
				else {
					if (temp_scat != 0) _M.PushBack(matrix_idx_local, group_matrix_idx, -volume * temp_scat);
				}
			}

			for (int idx = 1; idx <= 5; idx += 2) {
				if (adjacent_node_idx[idx] == INVALID) continue;
				_M.PushBack(matrix_idx_local, adjacent_node_idx[idx], -dTilda[idx]);
			}

			_chi_arr(ig, inode) = material._chi[ig];
		}
	}

	_M.RearrangeRow();
}

void Solver_t::Solve()
{
	Timer.TimeInit(Times::SOLVER_TIME);
	Control.Message("================================================================");
	Control.Message("    OUTER        keff            k diff.         psi diff.      ");
	Control.Message("================================================================");

	int icycle = 1;
	while (!_is_convergence) {

		UpdateRHS(icycle);

		LinearSolver(icycle);

		UpdateFissionSource(icycle);

		DetermineConvergence(icycle);

		icycle++;
	}

	GatherTotalFlux();

	GatherAssemblyFlux();

	Timer.TimeStop(Times::SOLVER_TIME);
}

void Solver_t::UpdateRHS(int icycle)
{
	double inv_keff = 1 / _keff;
	#pragma omp parallel for schedule(guided) collapse(2)
	for (int inode = 0; inode < _num_nodes; inode++) {
		for (int ig = 0; ig < _num_groups; ig++) {
			_RHS(ig, inode) = _chi_arr(ig, inode) * _psi(inode) * inv_keff;
		}
	}
}

void Solver_t::UpdateFissionSource(int icycle)
{	
	// Update fission source
	double volume = Geometry.GetUnitVolume();
	#pragma omp parallel for schedule(guided)
	for (int inode = 0; inode < _num_nodes; inode++) {
		int inode_global = inode + _node_offset;
		int imat = _mesh_materials(inode_global);
		Mat_t material = XS.GetMaterialInfo(imat);

		double sum = 0.0;
		for (int ig = 0; ig < _num_groups; ig++) {
			int idx = inode * _num_groups + ig + _row_offset;
			sum += material._nusigf[ig] * _phi(idx);
		}

		_psi(inode) = sum * volume;
	}

	// Update k-effective
	double temp_upper_dot = 0.0, temp_lower_dot = 0.0;
	for (int inode = 0; inode < _num_nodes; inode++) {
		temp_upper_dot += _psi(inode) * _psi(inode);
		temp_lower_dot += _psi(inode) * _prev_psi(inode);
	}

	double upper_dot = 0.0, lower_dot = 0.0;

	MPI_Reduce(&temp_upper_dot, &upper_dot, 1, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);
	MPI_Reduce(&temp_lower_dot, &lower_dot, 1, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);

	if (Control.GetMPIRank() == MASTER_PROCESS) _keff = _prev_keff * upper_dot / lower_dot;

	MPI_Bcast(&_keff, 1, MPI_DOUBLE, MASTER_PROCESS, MPI_COMM_WORLD);
}

void Solver_t::DetermineConvergence(int icycle)
{
	double keff_err = 1.0;
	double psi_err = 1.0, psi_sum = 1.0;

	double temp_err_sum = 0.0, temp_psi_sum = 0.0;

	for (int inode = 0; inode < _num_nodes; inode++) {
		temp_err_sum += (_psi(inode) - _prev_psi(inode)) * (_psi(inode) - _prev_psi(inode));
		temp_psi_sum += _psi(inode) * _psi(inode);
	}

	MPI_Reduce(&temp_err_sum, &psi_err, 1, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);
	MPI_Reduce(&temp_psi_sum, &psi_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);

	if ((Control.GetMPIRank() == MASTER_PROCESS)) {
		psi_sum = sqrt(psi_sum);
		psi_err = sqrt(psi_err) / psi_sum;

		keff_err = abs(_keff - _prev_keff) / abs(_keff);
		
		if ((psi_err < Control.GetPsiCritera()) && (keff_err < Control.GetEigenCritera())) _is_convergence = true;
	}

	MPI_Bcast(&_is_convergence, 1, MPI_C_BOOL, MASTER_PROCESS, MPI_COMM_WORLD);

	if (!_is_convergence) {
		_prev_keff = _keff;
		_prev_psi  = _psi;
	}

	stringstream s;
	s << "   " << icycle << "         " << _keff << "               " << keff_err << "               " << psi_err;
	Control.Message(s.str());
}

void Solver_t::GatherTotalFlux()
{
	int nz = Geometry.GetNumNodeZ();
	_total_phi.Create(_num_groups, nz * _unit_num_nodes);

	int* displs = new int[Control.GetMPISize()];
	int* counts = new int[Control.GetMPISize()];

	int* num_layers   = Geometry.GetNumLayers();
	int* start_layers = Geometry.GetStartLayers();

	for (int iproc = 0; iproc < Control.GetMPISize(); iproc++) {
		counts[iproc] = num_layers[iproc] * _unit_buffer_size;
		displs[iproc] = start_layers[iproc] * _unit_buffer_size;
	}

	_total_phi.GatherV(_phi, _row_offset, _num_rows, counts, displs);
}

void Solver_t::GatherAssemblyFlux()
{
	int num_horizon = 0, num_vertical = 0, num_axial = 0;
	int num_meshX = 0, num_meshY = 0, num_meshZ = 0;

	Geometry.GetBasicNum(num_horizon, num_vertical, num_axial);
	Geometry.GetNumMesh(num_meshX, num_meshY, num_meshZ);

	_assembly_phi.Create(_num_groups, num_horizon, num_vertical, num_axial); _assembly_phi = 0.0;

	int index = 0;
	for (int iz = 0; iz < num_axial; iz++) {
		for (int imeshZ = 0; imeshZ < num_meshZ; imeshZ++) {
			for (int iasyY = 0; iasyY < num_vertical; iasyY++) {
				for (int imeshY = 0; imeshY < num_meshY; imeshY++) {
					for (int iasyX = 0; iasyX < num_horizon; iasyX++) {
						if (Geometry.GetAssemblyIdx(iasyX, iasyY) == INVALID) continue;

						for (int imeshX = 0; imeshX < num_meshX; imeshX++) {
							for (int ig = 0; ig < _num_groups; ig++) 
								_assembly_phi(ig, iasyX, iasyY, iz) += _total_phi(ig, index);
							index++;
						}						
					}
				}
			}
		}
	}

	_assembly_phi /= num_meshX * num_meshY * num_meshZ;
}

void BiCGSTAB_t::Initialize()
{
	// flag
	_output_level = 1;
	_is_inner_convergence = false;

	// residual zero
	if (Control.GetMPIRank() != MASTER_PROCESS)
		_phi.IRecv(0, _unit_buffer_size, Control.GetMPIRank() - 1);
	if (Control.GetMPIRank() != Control.GetMPISize() - 1)
		_phi.IRecv(_row_offset + _num_rows, _unit_buffer_size, Control.GetMPIRank() + 1);

	if (Control.GetMPIRank() != MASTER_PROCESS)
		_phi.ISend(_row_offset, _unit_buffer_size, Control.GetMPIRank() - 1);
	if (Control.GetMPIRank() != Control.GetMPISize() - 1)
		_phi.ISend(_row_offset + _num_rows - _unit_buffer_size, _unit_buffer_size, Control.GetMPIRank() + 1);
	
	#pragma omp parallel for schedule(guided)
	for (int irow = 0; irow < _num_rows; irow++) {
		int num_cols  = _M.NumCol(irow);
		int start_col = _M.StartCol(irow);
		
		double sum = 0.0;
		for (int icol = 0; icol < num_cols; icol++) {
			int col = _M.ColIdx(start_col + icol);

			if ((col < _start_row) || (col > _end_row)) continue;
			sum += _phi(col - _start_row + _row_offset) * _M.Data(start_col + icol);
		}
		_temp_arr_buffer(irow) = sum;
	}

	_phi.WaitAll();
	
	#pragma omp parallel for schedule(guided)
	for (int irow = 0; irow < _num_rows; irow++) {
		int num_cols = _M.NumCol(irow);
		int start_col = _M.StartCol(irow);

		double sum = 0.0;
		for (int icol = 0; icol < num_cols; icol++) {
			int col = _M.ColIdx(start_col + icol);
			if ((col >= _start_row) && (col <= _end_row)) continue;

			sum += _phi(col - _start_row + _row_offset) * _M.Data(start_col + icol);
		}
		_temp_arr_buffer(irow) += sum;
	}

	_zero_residual = _RHS - _temp_arr_buffer;
	_residual      = _zero_residual;
	_prev_residual = _zero_residual;
	_hat_residual  = _zero_residual;

	// Calculate zero residual size
	double temp_dot_value = Dot(_zero_residual, _zero_residual);
	_zero_residual_size   = 0.0;

	MPI_Reduce(&temp_dot_value, &_zero_residual_size, 1, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);


	if (Control.GetMPIRank() == MASTER_PROCESS) _zero_residual_size = sqrt(_zero_residual_size);

	// Initialize scalar values
	_rho      = 1.0; _alpha  = 1.0; _omega = 1.0;
	_prev_rho = _rho; _prev_omega = _omega;

	// Initiliaze vector values
	_nu  = 0.0; _pvalue = 0.0;
	_prev_nu = _nu; _prev_pvalue = _pvalue;

	_hash   = 0.0; _svalue = 0.0; _tvalue = 0.0;
}

void BiCGSTAB_t::AllocateArray()
{
	_temp_arr_buffer.Create(_num_groups, _num_nodes);
	_zero_residual.Create(_num_groups, _num_nodes);
	_residual.Create(_num_groups, _num_nodes);
	_prev_residual.Create(_num_groups, _num_nodes);
	_hat_residual.Create(_num_groups, _num_nodes);

	_nu.Create(_num_groups, _num_nodes);      _pvalue.Create(_num_groups, _buffer_nodes_size);
	_prev_nu.Create(_num_groups, _num_nodes); _prev_pvalue.Create(_num_groups, _buffer_nodes_size);

	_hash.Create(_num_groups, _num_nodes);
	_svalue.Create(_num_groups, _buffer_nodes_size);
	_tvalue.Create(_num_groups, _num_nodes);
}

void BiCGSTAB_t::UpdatePvalue()
{
	// Rho calculation
	_temp_value_buffer = Dot(_hat_residual, _prev_residual);

	MPI_Reduce(&_temp_value_buffer, &_rho, 1, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);

	MPI_Bcast(&_rho, 1, MPI_DOUBLE, MASTER_PROCESS, MPI_COMM_WORLD);

	_beta = (_rho / _prev_rho) * (_alpha / _prev_omega);
	
	#pragma omp parallel for schedule(guided)
	for (int irow = 0; irow < _num_rows; irow++)
		_pvalue(irow + _row_offset) = _prev_residual(irow) + (_prev_pvalue(irow + _row_offset) - _prev_nu(irow) * _prev_omega) * _beta;
}

void BiCGSTAB_t::UpdateHash()
{
	// nu update

	if (Control.GetMPIRank() != MASTER_PROCESS)
		_pvalue.IRecv(0, _unit_buffer_size, Control.GetMPIRank() - 1);
	if (Control.GetMPIRank() != Control.GetMPISize() - 1)
		_pvalue.IRecv(_row_offset + _num_rows, _unit_buffer_size, Control.GetMPIRank() + 1);

	if (Control.GetMPIRank() != MASTER_PROCESS)
		_pvalue.ISend(_row_offset, _unit_buffer_size, Control.GetMPIRank() - 1);
	if (Control.GetMPIRank() != Control.GetMPISize() - 1)
		_pvalue.ISend(_row_offset + _num_rows - _unit_buffer_size, _unit_buffer_size, Control.GetMPIRank() + 1);
	
	#pragma omp parallel for schedule(guided)
	for (int irow = 0; irow < _num_rows; irow++) {
		int num_cols = _M.NumCol(irow);
		int start_col = _M.StartCol(irow);

		double sum = 0.0;
		for (int icol = 0; icol < num_cols; icol++) {
			int col = _M.ColIdx(start_col + icol);
			if ((col < _start_row) || (col > _end_row)) continue;

			sum += _pvalue(col - _start_row + _row_offset) * _M.Data(start_col + icol);
		}
		_nu(irow) = sum;
	}

	_pvalue.WaitAll();
	
	#pragma omp parallel for schedule(guided)
	for (int irow = 0; irow < _num_rows; irow++) {
		int num_cols = _M.NumCol(irow);
		int start_col = _M.StartCol(irow);

		double sum = 0.0;
		for (int icol = 0; icol < num_cols; icol++) {
			int col = _M.ColIdx(start_col + icol);
			if ((col >= _start_row) && (col <= _end_row)) continue;

			sum += _pvalue(col - _start_row + _row_offset) * _M.Data(start_col + icol);
		}
		_nu(irow) += sum;
	}

	// Update alpha
	_temp_value_buffer    = Dot(_hat_residual, _nu);

	double temp_dot_value = 0.0;
	
	MPI_Reduce(&_temp_value_buffer, &temp_dot_value, 1, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);

	MPI_Bcast(&temp_dot_value, 1, MPI_DOUBLE, MASTER_PROCESS, MPI_COMM_WORLD);

	_alpha = _rho / temp_dot_value;

	// Update Hash
	#pragma omp parallel for schedule(guided)
	for (int irow = 0; irow < _num_rows; irow++)
		_hash(irow) = _prev_phi(irow + _row_offset) + _alpha * _pvalue(irow + _row_offset);

}

void BiCGSTAB_t::UpdateResidual(int icycle)
{
	// Update s
	#pragma omp parallel for schedule(guided)
	for (int irow = 0; irow < _num_rows; irow++) {
		double multiple = _alpha * _nu(irow);
		_svalue(irow + _row_offset) = _prev_residual(irow) - multiple;
	}

	// Update t
	if (Control.GetMPIRank() != MASTER_PROCESS)
		_svalue.IRecv(0, _unit_buffer_size, Control.GetMPIRank() - 1);
	if (Control.GetMPIRank() != Control.GetMPISize() - 1)
		_svalue.IRecv(_row_offset + _num_rows, _unit_buffer_size, Control.GetMPIRank() + 1);

	if (Control.GetMPIRank() != MASTER_PROCESS)
		_svalue.ISend(_row_offset, _unit_buffer_size, Control.GetMPIRank() - 1);
	if (Control.GetMPIRank() != Control.GetMPISize() - 1)
		_svalue.ISend(_row_offset + _num_rows - _unit_buffer_size, _unit_buffer_size, Control.GetMPIRank() + 1);

	
	#pragma omp parallel for schedule(guided)
	for (int irow = 0; irow < _num_rows; irow++) {

		int num_cols  = _M.NumCol(irow);
		int start_col = _M.StartCol(irow);
		double sum = 0.0;
		for (int icol = 0; icol < num_cols; icol++) {
			int col = _M.ColIdx(start_col + icol);
			if (col < _start_row) continue;
			if (col > _end_row)   continue;

			sum += _svalue(col - _start_row + _row_offset) * _M.Data(start_col + icol);
		}

		_tvalue(irow) = sum;
	}

	_svalue.WaitAll();
	
	#pragma omp parallel for
	for (int irow = 0; irow < _num_rows; irow++) {

		int num_cols = _M.NumCol(irow);
		int start_col = _M.StartCol(irow);
		double sum = 0.0;
		for (int icol = 0; icol < num_cols; icol++) {
			int col = _M.ColIdx(start_col + icol);
			if ((col >= _start_row) && (col <= _end_row)) continue;

			sum += _svalue(col - _start_row + _row_offset) * _M.Data(start_col + icol);
		}

		_tvalue(irow) += sum;
	}

	// Update weight
	double temp_upper_dot = 0.0, temp_lower_dot = 0.0;

	for (int irow = 0; irow < _num_rows; irow++) {
		int row = irow + _row_offset;
		temp_upper_dot += _tvalue(irow) * _svalue(row);
		temp_lower_dot += _tvalue(irow) * _tvalue(irow);
	}
	
	double upper_dot = 0.0, lower_dot = 0.0;

	MPI_Reduce(&temp_upper_dot, &upper_dot, 1, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);
	MPI_Reduce(&temp_lower_dot, &lower_dot, 1, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);

	if (Control.GetMPIRank() == MASTER_PROCESS) _omega = upper_dot / lower_dot;

	MPI_Bcast(&_omega, 1, MPI_DOUBLE, MASTER_PROCESS, MPI_COMM_WORLD);

	// Update phi
	#pragma omp parallel 
	{
		#pragma omp for schedule(guided)
		for (int irow = 0; irow < _num_rows; irow++)
			_phi(irow + _row_offset) = _hash(irow) + _omega * _svalue(irow + _row_offset);

		// Update Residual
		#pragma	omp for schedule(guided)
		for (int irow = 0; irow < _num_rows; irow++)
			_residual(irow) = _svalue(irow + _row_offset) - _omega * _tvalue(irow);
	}

	// Determine convergence
	double temp_dot_buffer = Dot(_residual, _residual);
	double residual_size   = 0.0;

	MPI_Reduce(&temp_dot_buffer, &residual_size, 1, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);

	if (Control.GetMPIRank() == MASTER_PROCESS) {
		residual_size = sqrt(residual_size);
		double residual_error = residual_size / _zero_residual_size;
		if (residual_error < Control.GetInnerCriteria()) _is_inner_convergence = true;
	}

	MPI_Bcast(&_is_inner_convergence, 1, MPI_C_BOOL, MASTER_PROCESS, MPI_COMM_WORLD);

	_prev_phi = _phi;
	if (!_is_inner_convergence) {
		_prev_residual = _residual;
		_prev_rho      = _rho;
		_prev_omega    = _omega;
		_prev_nu       = _nu;
		_prev_pvalue   = _pvalue;
	}
}

void BiCGSTAB_t::LinearSolver(int icycle)
{
	Control.Message("      Solve linear system...");
	if (icycle == 1) AllocateArray();

	Initialize();
	 
	int iinner_cycle = 1;
	while (!_is_inner_convergence) {
	
		UpdatePvalue();
	
		UpdateHash();
	
		UpdateResidual(iinner_cycle);
	
		iinner_cycle++;
	}
}

/**
 * @brief Successive Over-relaxation method
 * @date 2022-06-03
 * parallel SOR algorithm using MPI and openmp
 * reference
 * (1) https://github.com/pranavladkat/SOR/blob/master/src/sor.hpp
 * (2) https://github.com/JasonPap/Parallel-SOR/blob/master/sor.c
 * (3) http://users.math.uoc.gr/~vagelis/Courses/Scientific_Computing/Papers/Parallel_SOR_Evans_84.pdf
 */

void SOR_t::Initialize()
{
	// flag
	_is_inner_convergence = false;

	// Initialize scalar values
	_w_optimal = 1.2;
	_spectral_radius = 0.0;
	_beta = 0.0;
	_diff_norm = 0.0;
	_is_inner_convergence = false;
	_max_inner_iterations = 128;
	_eps = 1e-8;
}

void SOR_t::AllocateArray()
{
	_matmul_buffer.Create(_num_groups, _num_nodes);
	_residual_buffer.Create(_num_groups, _num_nodes);
}

void SOR_t::CalculateWvalue()
{
	/** calculate optimal w value from given matrix
	 * @fixed me : calculate spectral radius of I-D^(-1)A
	 * reference : http://hpcf-files.umbc.edu/research/papers/YangGobbert2007SOR.pdf
	 */

	// rho : spectral radius of matrix I - D*(-1) * A
	// beta : rho(A), spectral radius of matrix I-D^(-1)A
	_w_optimal = 1.2;
}

double SOR_t::CalculateDiffNorm()
{
	double resid_per_procs = 0.0;
	int n_comp = _num_groups * _num_nodes;
	
	#pragma omp parallel
	{
		#pragma omp for
		for(int irow = 0; irow < _num_rows; irow++)
		{	
			resid_per_procs += (_phi(irow + _row_offset) - _prev_phi(irow + _row_offset)) * (_phi(irow + _row_offset) - _prev_phi(irow + _row_offset));
		}
	}

	resid_per_procs /= (double)n_comp;
	return sqrt(resid_per_procs);
}

void SOR_t::CheckIsConvergence()
{
	double resid_per_procs = CalculateDiffNorm();
	double resid = 0;

	MPI_Reduce(&resid_per_procs, &resid, 1, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);
	
	if (Control.GetMPIRank() == MASTER_PROCESS) {

		resid /= (double)Control.GetMPISize();
		_diff_norm = resid;

		if (_diff_norm < Control.GetInnerCriteria()){
			_is_inner_convergence = true;
		} 
	}

	MPI_Bcast(&_diff_norm, 1, MPI_DOUBLE, MASTER_PROCESS, MPI_COMM_WORLD);
	MPI_Bcast(&_is_inner_convergence, 1, MPI_C_BOOL, MASTER_PROCESS, MPI_COMM_WORLD);
	_prev_phi = _phi;
}

void SOR_t::UpdateVector()
{
	// Jacobi method with over-relaxation method(JOR)
	// phi(i)_{n+1} = (1-w) * phi(i)_n + w * (b - sum of a_ij * phi(j)_n) * 1 / a_ii

	if (Control.GetMPIRank() != MASTER_PROCESS)
	{
		_phi.IRecv(0, _unit_buffer_size, Control.GetMPIRank() - 1);
	}
		
	if (Control.GetMPIRank() != Control.GetMPISize() - 1)
	{
		_phi.IRecv(_row_offset + _num_rows, _unit_buffer_size, Control.GetMPIRank() + 1);
	}
		
	if (Control.GetMPIRank() != MASTER_PROCESS)
	{
		_phi.ISend(_row_offset, _unit_buffer_size, Control.GetMPIRank() - 1);
	}
		
	if (Control.GetMPIRank() != Control.GetMPISize() - 1)
	{
		_phi.ISend(_row_offset + _num_rows - _unit_buffer_size, _unit_buffer_size, Control.GetMPIRank() + 1);
	}

	int irow, icol;

	#pragma omp parallel 
	{
		#pragma omp for private(irow, icol)
		for (int irow = 0; irow < _num_rows; irow++)
		{
			int num_cols = _M.NumCol(irow);
			int start_col = _M.StartCol(irow);
			double sum = 0.0;

			for (int icol = 0; icol < num_cols; icol++)
			{
				int col = _M.ColIdx(start_col + icol);

				if ((col < _start_row) || (col > _end_row))
					continue;

				sum += _phi(col - _start_row + _row_offset) * _M.Data(start_col + icol);
			}

			_matmul_buffer(irow) = sum;
		}
	}

	_phi.WaitAll();

	#pragma omp parallel 
	{
		#pragma omp for private(irow, icol)
		for (int irow = 0; irow < _num_rows; irow++)
		{
			int num_cols = _M.NumCol(irow);
			int start_col = _M.StartCol(irow);
			double sum = 0.0;

			for (int icol = 0; icol < num_cols; icol++)
			{
				int col = _M.ColIdx(start_col + icol);

				if ((col >= _start_row) && (col <= _end_row))
					continue;
					
				sum += _phi(col - _start_row + _row_offset) * _M.Data(start_col + icol);
			}

			_matmul_buffer(irow) += sum;
		}
	}

	_residual_buffer = _RHS - _matmul_buffer;

	#pragma omp parallel
	{
		#pragma omp for private(irow)
		for(int irow = 0; irow < _num_rows; irow++)
		{	
			int start_col = _M.StartCol(irow);
			double a_ii = _M.Data(start_col + irow); // diagonal matrix a_ii

			_phi(irow + _row_offset) = _w_optimal * _residual_buffer(irow) / a_ii + (1 - _w_optimal) * _phi(irow + _row_offset);
		}
	}
}

void SOR_t::LinearSolver(int icycle)
{

	Control.Message("      Solve linear system...(SOR)");
	if (icycle == 1)
		AllocateArray();

	Initialize();
	CalculateWvalue();

	int iinner_cycle = 1;
	while (!_is_inner_convergence && iinner_cycle <= _max_inner_iterations)
	{
		UpdateVector();
		CheckIsConvergence();
		iinner_cycle++;
	}
}