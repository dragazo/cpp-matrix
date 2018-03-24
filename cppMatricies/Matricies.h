#ifndef MATRICIES_H
#define MATRICIES_H

#include <cstdlib>
#include <vector>
#include <exception>

struct MatrixSizeError : std::exception
{
	virtual const char *what() const override { return "Matrix had invalid size"; }
};

template<typename T>
class Matrix
{
private: // -- data -- //

	std::vector<T> data; // the elements in the array
	std::size_t r, c;    // number of rows / cols

public: // -- ctor / dtor / asgn -- //

	// creates an empty matrix
	Matrix() : r(0), c(0) {}

	// creates a rows x cols matrix
	Matrix(std::size_t rows, std::size_t cols) : r(rows), c(cols), data(rows * cols) {}

	// due to using std::vector for the data, default cpy ctor / dtor / asgn are sufficient

public: // -- operations -- //

	// gets the element at the specified row and col
	T &operator()(std::size_t row, std::size_t col) { return data[row * c + col]; }
	const T &operator()(std::size_t row, std::size_t col) const { return data[row * c + col]; }

	// gets the element at the specified row and col, but with bounds checking
	// throws std::out_of_range
	T &at(std::size_t row, std::size_t col)
	{
		if (row >= r || col >= c) throw std::out_of_range("Specified matrix element out of bounds");
		return data[row * c + col];
	}
	const T &at(std::size_t row, std::size_t col) const
	{
		if (row >= r || col >= c) throw std::out_of_range("Specified matrix element out of bounds");
		return data[row * c + col];
	}
};

// adds matrix b to matrix a
// throws MatrixSizeError if matricies are of different sizes
template<typename T>
static Matrix<T>& operator+=(Matrix<T> &a, const Matrix<T> &b)
{
	// ensure sizes are ok
	if (a.r != b.r || a.c != b.c) throw MatrixSizeError();

	for (std::size_t row = 0; row < r; ++row)
		for (std::size_t col = 0; col < c; ++col)
			a(row, col) += b(row, col);

	return a;
}
template<typename T>
static Matrix<T> operator+(const Matrix<T> &a, const Matrix<T> &b)
{
	Matrix<T> res = a;
	return res += b;
}

// subtracts matrix b from matrix a
// throws MatrixSizeError if matricies are of different sizes
template<typename T>
static Matrix<T>& operator-=(Matrix<T> &a, const Matrix<T> &b)
{
	// ensure sizes are ok
	if (a.r != b.r || a.c != b.c) throw MatrixSizeError();

	for (std::size_t row = 0; row < r; ++row)
		for (std::size_t col = 0; col < c; ++col)
			a(row, col) -= b(row, col);

	return a;
}
template<typename T>
static Matrix<T> operator-(const Matrix<T> &a, const Matrix<T> &b)
{
	Matrix<T> res = a;
	return res -= b;
}

#endif
