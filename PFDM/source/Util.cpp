#include "Util.h"

double Dot(Array_t<double> arr1, Array_t<double> arr2)
{
	double sum = 0.0;
	for (int i = 0; i < arr1.Size(); i++) sum += arr1(i) * arr2(i);
	return sum;
}

double Length(Array_t<double> arr1, Array_t<double> arr2)
{
	double sum = 0.0;
	for (int i = 0; i < arr1.Size(); i++) 
		sum += (arr1(i) - arr2(i)) * (arr1(i) - arr2(i));
	return sum;
}

template <typename T>
Array_t<T>& Array_t<T>::operator=(T val)
{
	for (int i = 0; i < _n; i++)
		_entry[i] = val;

	return *this;
}
template <typename T>
Array_t<T>& Array_t<T>::operator=(T* entry)
{
	for (int i = 0; i < _n; i++)
		_entry[i] = entry[i];

	return *this;
}

template <typename T>
Array_t<T>& Array_t<T>::operator=(const Array_t<T>& arr)
{
	if (!_is_allocate) Create(arr._nx, arr._ny, arr._nz, arr._nw);

	for (int i = 0; i < _n; i++)
		_entry[i] = arr._entry[i];

	return *this;
}

template <typename T>
Array_t<T>& Array_t<T>::operator*=(T val)
{
	for (int i = 0; i < _n; i++) _entry[i] = _entry[i] * val;

	return *this;
}

template <typename T>
Array_t<T>& Array_t<T>::operator/=(T val)
{
	for (int i = 0; i < _n; i++) _entry[i] = _entry[i] / val;

	return *this;
}

template <typename T>
Array_t<T> Array_t<T>::operator+(Array_t<T>& arr)
{
	Array_t<T> new_arr(arr._nx, arr._ny, arr._nz, arr._nw);
	
	for (int i = 0; i < _n; i++) new_arr._entry[i] = _entry[i] + arr._entry[i];

	return new_arr;
}

template <typename T>
Array_t<T> Array_t<T>::operator-(T val)
{
	Array_t<T> new_arr(this->_nx, this->_ny, this->_nz, this->_nw);
	
	for (int i = 0; i < _n; i++) new_arr._entry[i] = _entry[i] - val;

	return new_arr;
}

template <typename T>
Array_t<T> Array_t<T>::operator-(Array_t<T>& arr)
{
	Array_t<T> new_arr(arr._nx, arr._ny, arr._nz, arr._nw);
	
	for (int i = 0; i < _n; i++) new_arr._entry[i] = _entry[i] - arr._entry[i];

	return new_arr;
}

template <typename T>
Array_t<T> Array_t<T>::operator*(T val)
{
	Array_t<T> new_arr(this->_nx, this->_ny, this->_nz, this->_nw);
	
	for (int i = 0; i < _n; i++) new_arr._entry[i] = _entry[i] * val;

	return new_arr;
}

template <>
Array_t<double>& Array_t<double>::Normalize()
{
	double sum = 0.0;
	for (int i = 0; i < _n; i++) sum += _entry[i] * _entry[i];
	
	sum = sqrt(sum);

	*this /= sum;

	return *this;
}

template class Array_t<int>;
template class Array_t<double>;