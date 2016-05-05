//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <cstdlib>
#include <stdexcept>
#include <cassert>

template <class T>
class Matrix
{
public:
	Matrix(): _matrix(nullptr), _rows(0), _cols(0) {}
	Matrix(const Matrix& other):
		Matrix(other._rows, other._cols)
	{
		for (size_t i = 0; i < _rows; ++i)
			for (size_t j = 0; j < _cols; ++j)
				_matrix[i][j] = other._matrix[i][j];
	}
	Matrix(Matrix&& other):
		Matrix()
	{
		std::swap(_cols, other._cols);
		std::swap(_rows, other._rows);
		std::swap(_matrix, other._matrix);
	}
	Matrix& operator=(Matrix && other)
	{
		std::swap(_cols, other._cols);
		std::swap(_rows, other._rows);
		std::swap(_matrix, other._matrix);
		return *this;
	}
	Matrix& operator=(const Matrix& other)
	{
		Matrix temp(other);
		std::swap(_cols, temp._cols);
		std::swap(_rows, temp._rows);
		std::swap(_matrix, temp._matrix);
		return *this;
	}

	Matrix(size_t rows, size_t cols, T val = 0):
		_rows(rows), _cols(cols)
	{
		_matrix = new T*[rows];
		if (!rows)
			throw std::runtime_error("Zero matrix dimension");

		_matrix[0] = new T[rows * cols];
		 for (size_t i = 1; i < rows; ++i)
			_matrix[i] = _matrix[0] + i * cols;

		for (size_t i = 0; i < _rows; ++i)
			for (size_t j = 0; j < _cols; ++j)
				_matrix[i][j] = val;
	}
	~Matrix()
	{
		if (_matrix)
		{
			delete[] _matrix[0];
			delete[] _matrix;
		}
	}

	T& at(size_t row, size_t col) {return _matrix[row][col];}
	const T& at(size_t row, size_t col) const {return _matrix[row][col];}
	size_t nrows() const {return _rows;}
	size_t ncols() const {return _cols;}

private:
	T** _matrix;
	size_t _rows;
	size_t _cols;
};
