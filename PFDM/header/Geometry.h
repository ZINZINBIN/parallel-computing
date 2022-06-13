#pragma once
#include "Defines.h"
#include "Control.h"
#include "Util.h"
#include "Xsec.h"

class Geometry_t;

extern Geometry_t Geometry;

struct Assembly_t
{
	int _index;
	int _angle = 360;
	int* _axial_material;
};

class Geometry_t
{
private:
	int _n, _nx, _ny, _nz;
	int _base_nx, _base_ny, _base_nz;
	int _num_meshX, _num_meshY, _num_meshZ;   /// # of meshes per each pin
	int _num_mesh_per_assembly, _num_mesh_per_floor;
	int _num_assembly, _num_horizon, _num_vertical;
	int _rad_conf_angle;
	int *_start_layers, _start_layer, *_num_layers, _num_layer; /// MPI Parameters
	const int _num_adjacent_nodes = 6;
	bool _is_conf_cent = false;
	double _assembly_pitch;
	double _mesh_size[3], _surface_size[3];
	double _albedo[6], _boundary_diffusivity[6];
	double *_axial_node_size;
	double _hx, _hy, _hz;
	double _unit_volume;
	vector<Assembly_t> _assembly_types;
	Array_t<int> _rad_conf_mat, _rad_conf_idx;
	Array_t<int> _adjacent_nodes;
	Array_t<int> _mesh_materials;
	Array_t<double> _diffusivity;

	// Construct geometry configuration
	void SetGeometryParameters();
	void MappingGeometry();
	void MPISetting();

public:
	// Input setting functions
	void SetNumOfMesh(int num_meshX, int num_meshY, int num_meshZ) { 
		_num_meshX = num_meshX; _num_meshY = num_meshY; _num_meshZ = num_meshZ; 
	};
	void SetAssemblyPitch(double assembly_pitch)           { _assembly_pitch = assembly_pitch; };
	void SetNumOfAxialNode(int base_nz) { 
		_base_nz = base_nz;
		_axial_node_size = new double[_base_nz];
	}
	void SetAxialNode(int iz, double value)                { _axial_node_size[iz] = value; };
	void SetAlbedo(int i, double value)                    { _albedo[i] = value; };
	void SetAssemblyType(Assembly_t assembly)              { _assembly_types.push_back(assembly); };
	void SetRadConfAngle(int angle, bool is_cent)          { _rad_conf_angle = angle; _is_conf_cent = is_cent; };
	void SetRadConfSize(int num_assembly, int num_horizon, int num_vertical) { 
		_num_assembly = num_assembly; _num_horizon = num_horizon; _num_vertical = num_vertical;
		_rad_conf_mat.Create(_num_horizon + 2, _num_vertical + 2); _rad_conf_mat.Offset(-1, -1);
		_rad_conf_idx.Create(_num_horizon + 2, _num_vertical + 2); _rad_conf_idx.Offset(-1, -1);
	};
	void SetRadConfMat(Array_t<int>& rad_conf_mat)          { _rad_conf_mat = rad_conf_mat; };
	void SetRadConfIdx(Array_t<int>& rad_conf_idx)          { _rad_conf_idx = rad_conf_idx; };
	// Construct geometry configuration
	void ConstructGeometry();
	// Query functions
	int GetNumAdjacentNode()                           { return _num_adjacent_nodes; };
	int GetNumNodeZ()                                  { return _nz; };
	int GetNumMeshPerFloor()                           { return _num_mesh_per_floor; };
	int GetNumLayer()                                  { return _num_layer; };
	int GetStartLayer()                                { return _start_layer; };
	int GetAssemblyIdx(int iasyX, int iasyY)           { return _rad_conf_mat(iasyX, iasyY); };
	int* GetNumLayers()                                { return _num_layers; };
	int* GetStartLayers()                              { return _start_layers; };
	double GetAlbedo(int iadj)                         { return _albedo[iadj]; };
	double GetMeshSize(int iadj)                       { return _mesh_size[iadj]; };
	double GetUnitVolume()                             { return _unit_volume; };
	/// Transfer info to other classes
	void GetAdjacentNode(Array_t<int>& adjacent_nodes) { adjacent_nodes = _adjacent_nodes; };
	void GetMeshMaterial(Array_t<int>& mesh_materials) { mesh_materials = _mesh_materials; };
	void GetSurfaceSize(double* surface_size)           { for (int i = 0; i < 3; i++) surface_size[i] = _surface_size[i]; };
	void GetDiffusivity(Array_t<double>& diffusivity) { diffusivity = _diffusivity; };
	void GetBoundaryDiffusivity(double* boundary_diffusivity) {
		for (int i = 0; i < 6; i++) boundary_diffusivity[i] = _boundary_diffusivity[i];
	};
	void GetNumMesh(int& num_meshX, int& num_meshY, int& num_meshZ) {
		num_meshX = _num_meshX; num_meshY = _num_meshY; num_meshZ = _num_meshZ;
	};
	void GetBasicNum(int& num_horizon, int& num_vertical, int& num_axial) {
		num_horizon = _num_horizon; num_vertical = _num_vertical; num_axial = _base_nz;
	};
};
