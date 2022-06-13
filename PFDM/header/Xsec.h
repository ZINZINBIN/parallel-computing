#pragma once
#include "Defines.h"
#include "Util.h"

/// Xsec.h
class XS_t;

extern XS_t XS;

struct Mat_t
{
	int _index;
	int _ng;
	double *_D;
	double *_sigr;
	double *_nusigf;
	double *_chi;
	Array_t<double> _scat;
	void SetNumOfGroup(int ng) {
		_ng     = ng;
		_D      = new double[ng];
		_sigr   = new double[ng];
		_nusigf = new double[ng];
		_chi    = new double[ng];
		_scat.Create(ng, ng);
	};
	void Destroy() { delete[] _D, _sigr, _nusigf, _chi, _scat; };
};

class XS_t
{
private:
	int _lib_type;
	int _ng;
	int _num_material;
	vector<Mat_t> _mat_types;
	
public:
	// Setting functions
	void SetLibType(int lib_type)        { _lib_type = lib_type; };
	void SetNumGroup(int ng)             { _ng = ng; };
	void SetMaterialType(Mat_t mat_type) { _mat_types.push_back(mat_type); };

	// Query functions
	int GetNumOfGruop()                  { return _ng; };
	Mat_t GetMaterialInfo(int idx)       { return _mat_types[idx]; };
};