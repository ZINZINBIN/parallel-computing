#pragma once
#include "Defines.h"

double Dot(Array_t<double> arr1, Array_t<double> arr2);
double Length(Array_t<double> arr1, Array_t<double> arr2);

template <typename T>
class Array_t
{
private:
	T* _entry;
	int _n = 1, _nx = 1, _ny = 1, _nz = 1, _nw = 1;
	int _nxy = 1, _nxyz = 1;
	int _dx = 0, _dy = 0, _dz = 0, _dw = 0;
	int _recv_count;
	bool _is_allocate = false;
	vector<MPI_Request> _requests;

public:
	~Array_t() { Destroy(); };
	Array_t() {};
	Array_t(const Array_t<T>& arr) {
		this->operator=(arr);
		this->Offset(arr._dx, arr._dy, arr._dz, arr._dw);
	};
	Array_t(int nx, int ny=1, int nz=1, int nw=1)         { Create(nx, ny, nz, nw); };
	Array_t(Array_t&& arr);
	void Destroy();
	void Create(int nx, int ny=1, int nz=1, int nw=1);
	void Offset(int dx, int dy=0, int dz=0, int dw=0);
	int  Size()                                            { return _n; };
	T* Pointer()                                           { return _entry; };
	/// MPI functions
	void ISend(int start, int count, int trg_rank);
	void IRecv(int start, int count, int src_rank);
	void WaitAll();
	void GatherV(Array_t<T>& send_arr, int offset, int count, int* counts, int* displs);
	/// Operators
	inline T& operator() (int ix)                         { return _entry[ix - _dx]; };
	inline T& operator() (int ix, int iy)                 { return _entry[ix - _dx + _nx * (iy - _dy)]; };
	inline T& operator() (int ix, int iy, int iz)         { return _entry[ix - _dx + _nx * (iy - _dy) + _nxy * (iz - _dz)]; };
	inline T& operator() (int ix, int iy, int iz, int iw) { return _entry[ix - _dx + _nx * (iy - _dy) + _nxy * (iz - _dz) + _nxyz * (iw - _dx)]; };
	inline T& operator[] (int ix)                         { return _entry[ix - _dx]; };
	Array_t<T>& operator=(T val);
	Array_t<T>& operator=(T* entry);
	Array_t<T>& operator=(const Array_t<T>& arr);
	Array_t<T>& operator*=(T val);
	Array_t<T>& operator/=(T val);
	Array_t<T> operator+(Array_t<T>& arr);
	Array_t<T> operator-(T val);
	Array_t<T> operator-(Array_t<T>& arr);
	Array_t<T> operator*(T val);
	Array_t<T>& Normalize();
};

template <typename T>
void Array_t<T>::Destroy()
{
	if (!_is_allocate) return;

	if (_entry != nullptr) delete[] _entry;

	_is_allocate = false;
}

template <typename T>
void Array_t<T>::Create(int nx, int ny, int nz, int nw)
{
	if (_is_allocate) Destroy();

	_nx   = nx;
	_ny   = ny;
	_nz   = nz;
	_nw   = nw;
	_nxy  = _nx * _ny;
	_nxyz = _nxy * _nz;
	_n    = _nxyz * _nw;
	
	_is_allocate = true;
	_entry = new T[_n];
}

template <typename T>
Array_t<T>::Array_t(Array_t&& arr)
{
	this->_nx   = arr._nx;
	this->_ny   = arr._ny;
	this->_nz   = arr._nz;
	this->_nw   = arr._nw;
	this->_nxy  = arr._nxy;
	this->_nxyz = arr._nxyz;
	
	if (_is_allocate) delete[] this->_entry;

	this->_entry = arr._entry;

	arr._entry = nullptr;
}

template <typename T>
void Array_t<T>::Offset(int dx, int dy, int dz, int dw)
{
	_dx = dx; _dy = dy; _dz = dz; _dw = dw;
}

template <typename T>
void Array_t<T>::ISend(int start, int count, int src_rank)
{
	int tag = 0;
	MPI_Request request;

	MPI_Isend(_entry + start, count, MPI_DOUBLE, src_rank, tag, MPI_COMM_WORLD, &request);
	
	_requests.push_back(request);
}

template <typename T>
void Array_t<T>::IRecv(int start, int count, int src_rank)
{
	int tag = 0;
	MPI_Request request;

	MPI_Irecv(_entry + start, count, MPI_DOUBLE, src_rank, tag, MPI_COMM_WORLD, &request);

	_requests.push_back(request);
}

template <typename T>
void Array_t<T>::WaitAll()
{
	MPI_Waitall(_requests.size(), _requests.data(), MPI_STATUSES_IGNORE);

	_requests.clear();
}

template <typename T>
void Array_t<T>::GatherV(Array_t<T>& send_arr, int offset, int count, int* counts, int* displs)
{
	MPI_Gatherv(send_arr.Pointer() + offset, count, MPI_DOUBLE, _entry, counts, displs, MPI_DOUBLE, MASTER_PROCESS, MPI_COMM_WORLD);
}

template <typename T>
class CSR_t
{
private:
	int _num_row = 0;
	int *_row;
	bool _is_row_allocate = false;
	bool _is_element_clear = true;
	vector<int> _col;
	vector<T>   _data;

public:
	void SetRowSize(int num_row);
	void Destroy();
	void ClearElements();
	void PushBack(int irow, int icol, T val);
	void RearrangeRow();
	int Size()              { return _row[_num_row]; };
	int StartCol(int irow)  { return _row[irow]; };
	int NumCol(int irow)    { return _row[irow + 1] - _row[irow]; };

	int ColIdx(int idx)     { return _col[idx]; };
	T&  Data(int idx)       { return _data[idx]; };
};

template <typename T>
void CSR_t<T>::SetRowSize(int num_row)
{
	Destroy();

	_num_row = num_row;
	_row     = new int[num_row + 1];

	memset(_row, 0, (_num_row + 1) * sizeof(int));
}

template <typename T>
void CSR_t<T>::Destroy()
{
	if (!_is_row_allocate) return;

	ClearElements();

	_is_row_allocate = false;
	_num_row         = 0;

	delete[] _row;
}

template <typename T>
void CSR_t<T>::ClearElements()
{
	if (_is_element_clear) return;

	_is_element_clear = true;

	_col.clear();  vector<int>().swap(_col);
	_data.clear(); vector<T>().swap(_data);

	memcpy(_row, 0, (_num_row + 1) * sizeof(int));
}

template <typename T>
void CSR_t<T>::PushBack(int irow, int icol, T val)
{
	_row[irow + 1]++;

	_col.push_back(icol);
	_data.push_back(val);
}

template <typename T>
void CSR_t<T>::RearrangeRow()
{
	for (int irow = 1; irow < _num_row + 1; irow++) _row[irow] += _row[irow - 1];
}