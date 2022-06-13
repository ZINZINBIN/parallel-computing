#include "Geometry.h"

Geometry_t Geometry;

void Geometry_t::ConstructGeometry()
{
	Timer.TimeInit(Times::GEOMETRY_TIME);
	
	MappingGeometry();

	SetGeometryParameters();

	MPISetting();

	Timer.TimeStop(Times::GEOMETRY_TIME);
}

void Geometry_t::MappingGeometry()
{
	Control.Message("Mapping geometry cells...");
	// Total number of mesh in core
	_num_mesh_per_assembly = _num_meshX * _num_meshY;
	_num_mesh_per_floor    = _num_meshX * _num_meshY * _num_assembly;
	_nz                    = _num_meshZ * _base_nz;
	_n                     = _num_mesh_per_floor * _nz;

	// Set adjacent assembly mapping
	int *num_horizons = new int[_num_vertical];
	for (int iy = 0; iy < _num_vertical; iy++) {
		int num_horizon = 0;
		for (int ix = 0; ix < _num_horizon; ix++) {
			if (_rad_conf_idx(ix, iy) != INVALID) num_horizon++;
		}
		num_horizons[iy] = num_horizon;
	}
	
	// Set adjacent mapping array and Allocate materials in each mesh
	_adjacent_nodes.Create(_num_adjacent_nodes, _n);
	_mesh_materials.Create(_n);
	_adjacent_nodes = INVALID; _mesh_materials = INVALID;

	int index = 0;
	for (int iz = 0; iz < _base_nz; iz++) {
		for (int imeshZ = 0; imeshZ < _num_meshZ; imeshZ++) {
			for (int iasyY = 0; iasyY < _num_vertical; iasyY++) {
				for (int imeshY = 0; imeshY < _num_meshY; imeshY++) {
					for (int iasyX = 0; iasyX < _num_horizon; iasyX++) {
						if (_rad_conf_idx(iasyX, iasyY) == INVALID) continue;
						for (int imeshX = 0; imeshX < _num_meshX; imeshX++) {
							// Allocate material
							int iasy = _rad_conf_mat(iasyX, iasyY);
							int imat = _assembly_types[iasy]._axial_material[iz];
							_mesh_materials(index) = imat;

							// Adjacent nodes
							/// West
							if ((imeshX == 0) && (_rad_conf_idx(iasyX - 1, iasyY) == INVALID))
								_adjacent_nodes(WEST, index) = INVALID;
							else
								_adjacent_nodes(WEST, index) = index - 1;

							/// East
							if ((imeshX == _num_meshX - 1) && (_rad_conf_idx(iasyX + 1, iasyY) == INVALID))
								_adjacent_nodes(EAST, index) = INVALID;
							else
								_adjacent_nodes(EAST, index) = index + 1;

							/// North
							if (imeshY == 0) {
								if (_rad_conf_idx(iasyX, iasyY - 1) == INVALID)
									_adjacent_nodes(NORTH, index) = INVALID;
								else
									_adjacent_nodes(NORTH, index) = index - _num_meshX * (_rad_conf_idx(iasyX, iasyY) - _rad_conf_idx(iasyX, iasyY - 1));
							}
							else    _adjacent_nodes(NORTH, index) = index - _num_meshX * num_horizons[iasyY];

							/// South
							if (imeshY == _num_meshY - 1) {
								if (_rad_conf_idx(iasyX, iasyY + 1) == INVALID)
									_adjacent_nodes(SOUTH, index) = INVALID;
								else
									_adjacent_nodes(SOUTH, index) = index + _num_meshX * (_rad_conf_idx(iasyX, iasyY + 1) - _rad_conf_idx(iasyX, iasyY));
							}
							else    _adjacent_nodes(SOUTH, index) = index + _num_meshX * num_horizons[iasyY];

							/// Lower
							if ((iz == 0) && (imeshZ == 0)) 
								_adjacent_nodes(LOWER, index) = INVALID;
							else 
								_adjacent_nodes(LOWER, index) = index - _num_mesh_per_floor;

							/// Upper
							if ((iz == _base_nz - 1) && (imeshZ == _num_meshZ - 1)) 
								_adjacent_nodes(UPPER, index) = INVALID;
							else 
								_adjacent_nodes(UPPER, index) = index + _num_mesh_per_floor;
							index++;
						}
					}
				}
			}
		}
	}

	delete[] num_horizons;
}

void Geometry_t::SetGeometryParameters()
{
	_hx = _assembly_pitch / _num_meshX;
	_hy = _assembly_pitch / _num_meshY;
	_hz = _axial_node_size[0] / _num_meshZ;

	_mesh_size[0] = _hx;
	_mesh_size[1] = _hy;
	_mesh_size[2] = _hz;

	_surface_size[0] = _hy * _hz;
	_surface_size[1] = _hx * _hz;
	_surface_size[2] = _hy * _hx;

	_unit_volume = _hx * _hy * _hz;

	for (int idim = 0; idim < 3; idim++) {

		for (int idir = 0; idir < 2; idir++)
			_boundary_diffusivity[2 * idim + idir] = _albedo[2 * idim + idir] / 2.0;
	}

	// Diffusivity
	_diffusivity.Create(3, XS.GetNumOfGruop(), _n);

	double mesh = 0.0;
	double diffusion_coefficient = 0.0;
	for (int inode = 0; inode < _n; inode++) {
		int imat = _mesh_materials(inode);
		Mat_t material = XS.GetMaterialInfo(imat);

		for (int ig = 0; ig < XS.GetNumOfGruop(); ig++) {
			for (int idir = 0; idir < 3; idir++) {
				mesh = _mesh_size[idir];

				_diffusivity(idir, ig, inode) = material._D[ig] / mesh;
			}
		}
	}
}

void Geometry_t::MPISetting()
{
	Control.Message("Setting MPI environment...");
	_start_layers = new int[Control.GetMPISize()];
	_num_layers   = new int[Control.GetMPISize()];

	if (Control.GetMPIRank() == MASTER_PROCESS) {
		int share    = (int) _nz / Control.GetMPISize();
		int residual = (int) _nz % Control.GetMPISize();

		int sum = 0;
		for (int i = 0; i < Control.GetMPISize(); i++) {
			_num_layers[i] = share;
			if (i < residual) _num_layers[i]++;

			_start_layers[i] = sum;
			sum += _num_layers[i];
		}
	}

	MPI_Bcast(_start_layers, Control.GetMPISize(), MPI_INTEGER, MASTER_PROCESS, MPI_COMM_WORLD);
	MPI_Bcast(_num_layers, Control.GetMPISize(), MPI_INTEGER, MASTER_PROCESS, MPI_COMM_WORLD);

	_start_layer = _start_layers[Control.GetMPIRank()];
	_num_layer   = _num_layers[Control.GetMPIRank()];
}