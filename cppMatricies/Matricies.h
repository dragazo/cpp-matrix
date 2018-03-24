#ifndef MATRICIES_H
#define MATRICIES_H

#include <cstdlib>
#include <vector>
#include <exception>
#include <utility>

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
	Matrix(std::size_t rows, std::size_t cols) : r(rows), c(cols), data(rows * cols)
	{
		// if either dimension was zero, both are zero
		if (rows == 0 || cols == 0) r = c = 0;
	}

	// due to using std::vector for the data, default cpy ctor / dtor / asgn are sufficient

public: // -- operations -- //

	// gets the number of rows/cols in the matrix
	std::size_t rows() const { return r; }
	std::size_t cols() const { return c; }

	// returns the number of elements in the matrix
	std::size_t size() const { return data.size(); }
	// returns the current capacity of the matrix
	std::size_t capacity() const { return data.capacity(); }

	// returns true if the matrix is square
	bool square() const { return r == c; }
	// returns true if the matrix is empty
	bool empty() const { return r == 0; }

	// resizes the matrix
	void resize(std::size_t rows, std::size_t cols)
	{
		// if either dimension was zero, both are zero
		if (rows == 0 || cols == 0)
		{
			r = c = 0;
			data.resize(0);
		}
		// otherwise resize as usual
		else
		{
			r = rows;
			c = cols;
			data.resize(rows * cols);
		}
	}
	// requests the matrix to set aside space for the specified number of elements
	void reserve(std::size_t count) { data.reserve(count); }

	// requests the matrix to remove unused allocated space (request may be declined)
	void shrink_to_fit() { data.shrink_to_fit(); }

	// gets the element at the specified row and col
	T &operator()(std::size_t row, std::size_t col) { return data[row * c + col]; }
	const T &operator()(std::size_t row, std::size_t col) const { return data[row * c + col]; }

	// gets the element at the specified row and col, but with bounds checking
	// throws std::out_of_range if element is out of bounds
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
static Matrix<T> &operator+=(Matrix<T> &a, const Matrix<T> &b)
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
static Matrix<T> &operator-=(Matrix<T> &a, const Matrix<T> &b)
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

// multiplies the matrix by a scalar
template<typename T>
static Matrix<T> &operator*=(Matrix<T> &m, const T &f)
{
	for (std::size_t row = 0; row < r; ++row)
		for (std::size_t col = 0; col < c; ++col)
			a(row, col) *= f;
}
template<typename T>
static Matrix<T> operator*(const Matrix<T> &m, const T &f)
{
	Matrix<T> res = m;
	return res *= f;
}
template<typename T>
static Matrix<T> operator*(const T &f, const Matrix<T> &m)
{
	Matrix<T> res = m;
	return res *= f;
}

// multiplies matrix lhs by matrix rhs
// throws MatrixSizeError if matricies are incompatible
template<typename T>
static Matrix<T> operator*(const Matrix<T> &lhs, const Matrix<T> &rhs)
{
	// matricies must be compatible
	if (lhs.cols() != rhs.rows()) throw MatrixSizeError();

	Matrix<T> res(lhs.rows(), rhs.cols()); // allocate the result
	T dot; // the destination of the dot product (potentially-large)

	for(std::size_t row = 0; row < res.rows(); ++row)
		for (std::size_t col = 0; col < res.cols(); ++col)
		{
			// initialize dot to first term (as T is not guaranteed to have a ctor with a default of "zero")
			dot = lhs(row, 0) * rhs(0, col);

			// add all the other terms
			for (std::size_t i = 1; i < lhs.cols(); ++i)
				dot += lhs(row, i) * rhs(i, col);

			// assign as result entry
			res(row, col) = std::move(dot);
		}

	return res;
}
template<typename T>
static Matrix<T> operator*=(Matrix<T> &lhs, const Matrix<T> &rhs)
{
	return *this = lhs * rhs;
}

#endif
