#ifndef MATRICIES_H
#define MATRICIES_H

#include <cstdlib>
#include <vector>
#include <exception>
#include <utility>
#include <iostream>

struct MatrixSizeError : std::exception
{
private:
	const char *msg;

public:
	MatrixSizeError() : msg("") {}
	MatrixSizeError(const char *_msg) : msg(_msg) {}

	virtual const char *what() const override { return msg; }
};

template<typename T>
class Matrix
{
public: // -- enums / etc -- //

	// the result of a row reduction operation
	/*enum RRResult
	{
		None = 0, REF = 1, RREF = 3
	};*/

private: // -- data -- //

	std::vector<T> data; // the elements in the array
	std::size_t r, c;    // number of rows / cols

private: // -- helpers -- //



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

public: // -- utilities -- //

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
	T operator()(std::size_t row, std::size_t col) const { return data[row * c + col]; }

	// gets the element at the specified row and col, but with bounds checking
	// throws std::out_of_range if element is out of bounds
	T &at(std::size_t row, std::size_t col)
	{
		if (row >= r || col >= c) throw std::out_of_range("Specified matrix element out of bounds");
		return data[row * c + col];
	}
	T at(std::size_t row, std::size_t col) const
	{
		if (row >= r || col >= c) throw std::out_of_range("Specified matrix element out of bounds");
		return data[row * c + col];
	}

	// creates an nxn identity matrix
	static Matrix identity(std::size_t n)
	{
		// allocate the result
		Matrix res(n, n);

		// fill with identity elements
		for (std::size_t row = 0; row < n; ++row)
			for (std::size_t col = 0; col < n; ++col)
				res(row, col) = row == col ? 1 : 0;

		return res;
	}

public: // -- elementary row operations -- //

	// swaps rows a and b
	void swapRows(std::size_t a, std::size_t b)
	{
		using std::swap; // ADL idiom

		// early exit if swapping to same destination
		if (a == b) return;

		for (std::size_t i = 0; i < c; ++i)
			swap((*this)(a, i), (*this)(b, i));
	}

	// multiplies the elements in a row by a scalar
	void multRow(std::size_t row, T f)
	{
		for (std::size_t i = 0; i < c; ++i) (*this)(row, i) *= f;
	}
	// divides the elements in a row by a scalar
	void divRow(std::size_t row, T f)
	{
		for (std::size_t i = 0; i < c; ++i) (*this)(row, i) /= f;
	}

	// adds row <from> to row <to>
	void addRow(std::size_t to, std::size_t from)
	{
		for (std::size_t i = 0; i < c; ++i) (*this)(to, i) += (*this)(from, i);
	}
	// subtracts row <from> from row <to>
	void subRow(std::size_t to, std::size_t from)
	{
		for (std::size_t i = 0; i < c; ++i) (*this)(to, i) -= (*this)(from, i);
	}

	// adds a multiple of row <from> to row <to>
	void addRowMult(std::size_t to, std::size_t from, T f)
	{
		for (std::size_t i = 0; i < c; ++i) (*this)(to, i) += (*this)(from, i) * f;
	}
	// subtracts a multiple of row <from> to row <to>
	void subRowMult(std::size_t to, std::size_t from, T f)
	{
		for (std::size_t i = 0; i < c; ++i) (*this)(to, i) -= (*this)(from, i) * f;
	}

public: // -- operations -- //

	// finds the determinant of the matrix
	// throws MatrixSizeError if cannot find determinant or if matrix is empty
	T det() const
	{
		if (r != c) throw MatrixSizeError("Only square matricies have a determinant");
		if (r == 0) throw MatrixSizeError("Cannot take the determinant of an empty matrix");

		switch (r)
		{
		case 1: return (*this)(0, 0);
		case 2: return (*this)(0, 0) * (*this)(1, 1) - (*this)(1, 0)*(*this)(0, 1);

		default:
			T sum = 0; // initialize sum to zero

			// add cofactors along first row
			for (std::size_t col = 0; col < c; ++col)
				sum += cofactor(0, col);

			return sum;
		}
	}

	// returns the <row><col> minor (i.e. the result of removing row <row> and col <col>)
	// throws MatrixSizeError if matrix is empty
	// also throws by det()
	Matrix minor(std::size_t row, std::size_t col) const
	{
		if (r == 0) throw MatrixSizeError("Cannot take a minor from an empty matrix");

		// allocate the result
		Matrix res(r - 1, c - 1);

		// fill its entries
		for (std::size_t _row = 0; _row < row; ++_row)
		{
			for (std::size_t _col = 0; _col < col; ++_col)
				res(_row, _col) = (*this)(_row, _col);
			for (std::size_t _col = col + 1; _col < c; ++_col)
				res(_row, _col - 1) = (*this)(_row, _col);
		}
		for (std::size_t _row = row + 1; _row < r; ++_row)
		{
			for (std::size_t _col = 0; _col < col; ++_col)
				res(_row - 1, _col) = (*this)(_row, _col);
			for (std::size_t _col = col + 1; _col < c; ++_col)
				res(_row - 1, _col - 1) = (*this)(_row, _col);
		}

		return res;
	}
	// returns the <row><col> cofactor (i.e. (-1)^(row+col) * (row, col) * minor(row,col).det() )
	// throws by minor()
	T cofactor(std::size_t row, std::size_t col) const
	{
		return ((row + col) & 1 ? -1 : 1) * (*this)(row, col) * minor(row, col).det();
	}
	// returns the adjugate matrix (i.e. matrix of cofactors)
	// throws by cofactor()
	Matrix adjugate() const
	{
		// allocate the result
		Matrix res(r, c);

		// fill with cofactors
		for (std::size_t row = 0; row < r; ++row)
			for (std::size_t col = 0; col < c; ++col)
				res(row, col) = cofactor(row, col);

		return res;
	}

	// puts the matrix into row echelon form. returns the rank of the matrix
	std::size_t REF()
	{
		// for each row
		for (std::size_t row = 0, col; row < r; ++row)
		{
			// find the leading entry
			for (col = row; col < c; ++col)
			{
				// move all the zero-entries here and down in this col to the bottom
				for (std::size_t top = row, bottom = r - 1; ;)
				{
					// wind top to next zero entry
					for (; top < bottom && (*this)(top, col) != 0; ++top);
					// wind bottom to next nonzero entry
					for (; top < bottom && (*this)(bottom, col) == 0; --bottom);

					// if they're still valid, perform the swap
					if (top < bottom) swapRows(top, bottom);
					// otherwise, we're done sorting
					else break;
				}

				// if after sorting this element is nonzero, it is the leading entry
				if ((*this)(row, col) != 0) break;
			}
			// if we didn't find a leading entry, the rest of the matrix is zeroes (so we're done)
			if(col == c) return row;

			// make this row's leading entry a 1 via row division
			divRow(row, (*this)(row, col));

			// use this row to eliminate the lower rows
			for (std::size_t j = row + 1; j < r && (*this)(j, col) != 0; ++j)
				subRowMult(j, row, (*this)(j, col));
		}

		// successfully converted to REF
		return r;
	}
	// puts the matrix into reduced row echelon form. returns the rank of the matrix
	std::size_t RREF()
	{
		// for each row
		for (std::size_t row = 0, col; row < r; ++row)
		{
			// find the leading entry
			for (col = row; col < c; ++col)
			{
				// move all the zero-entries here and down in this col to the bottom
				for (std::size_t top = row, bottom = r - 1; ;)
				{
					// wind top to next zero entry
					for (; top < bottom && (*this)(top, col) != 0; ++top);
					// wind bottom to next nonzero entry
					for (; top < bottom && (*this)(bottom, col) == 0; --bottom);

					// if they're still valid, perform the swap
					if (top < bottom) swapRows(top, bottom);
					// otherwise, we're done sorting
					else break;
				}

				// if after sorting this element is nonzero, it is the leading entry
				if ((*this)(row, col) != 0) break;
			}
			// if we didn't find a leading entry, the rest of the matrix is zeroes (so we're done)
			if (col == c) return row;

			// make this row's leading entry a 1 via row division
			divRow(row, (*this)(row, col));

			// use this row to eliminate the higher rows
			for (std::size_t j = 0; j < row; ++j)
				if ((*this)(j, col) != 0) subRowMult(j, row, (*this)(j, col));

			// use this row to eliminate the lower rows
			for (std::size_t j = row + 1; j < r && (*this)(j, col) != 0; ++j)
				subRowMult(j, row, (*this)(j, col));
		}

		// successfully converted to REF
		return r;
	}

	// attempts to find the inverse of the matrix. returns true if successful and stores inverse in <dest>
	// throws MatrixSizeError if matrix is not square or is empty
	bool inverse(Matrix &dest) const
	{
		if (r != c) throw MatrixSizeError("Only square matricies are invertible");
		if (r == 0) throw MatrixSizeError("Cannot take the inverse of an empty matrix");

		Matrix util(r, 2 * c); // make a wide matrix

		// fill left side with current matrix and right side with identity matrix
		for (std::size_t row = 0; row < r; ++row)
			for (std::size_t col = 0; col < c; ++col)
			{
				util(row, col) = (*this)(row, col);
				util(row, col + c) = row == col ? 1 : 0;
			}

		// if putting util matrix into rref has full rank, inversion was successful
		if (util.RREF() == r)
		{
			// allocate dest
			dest.resize(r, c);

			// copy inverse (right side) to dest
			for (std::size_t row = 0; row < r; ++row)
				for (std::size_t col = 0; col < c; ++col)
					dest(row, col) = util(row, col + c);

			return true;
		}
		else return false;
	}
};

// -- io -- //

template<typename T>
std::ostream &operator<<(std::ostream &ostr, const Matrix<T> &m)
{
	for (std::size_t row = 0; row < m.rows(); ++row)
	{
		for (std::size_t col = 0; col < m.cols(); ++col)
			ostr << std::setw(10) << m(row, col) << ' ';

		ostr << '\n';
	}

	return ostr;
}

// -- operator definitions -- //

// adds matrix b to matrix a
// throws MatrixSizeError if matricies are of different sizes
template<typename T>
static Matrix<T> &operator+=(Matrix<T> &a, const Matrix<T> &b)
{
	// ensure sizes are ok
	if (a.r != b.r || a.c != b.c) throw MatrixSizeError("Matrix addition requires the matricies be of same size");

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
	if (a.r != b.r || a.c != b.c) throw MatrixSizeError("Matrix subtraction requires the matricies be of same size");

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
	if (lhs.cols() != rhs.rows()) throw MatrixSizeError("Matrix multiplication requires lhs #cols equal rhs #rows");

	Matrix<T> res(lhs.rows(), rhs.cols()); // allocate the result
	T dot; // the destination of the dot product

	for(std::size_t row = 0; row < res.rows(); ++row)
		for (std::size_t col = 0; col < res.cols(); ++col)
		{
			// initialize dot to zero
			dot = 0;

			// compute the dot product
			for (std::size_t i = 0; i < lhs.cols(); ++i)
				dot += lhs(row, i) * rhs(i, col);

			// assign as result entry
			res(row, col) = dot;
		}

	return res;
}
template<typename T>
static Matrix<T> operator*=(Matrix<T> &lhs, const Matrix<T> &rhs)
{
	return *this = lhs * rhs;
}

#endif
