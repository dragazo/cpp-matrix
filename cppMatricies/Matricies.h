#ifndef MATRICIES_H
#define MATRICIES_H

#include <cstdlib>
#include <vector>
#include <exception>
#include <utility>
#include <iostream>
#include <type_traits>

// macros for accessing <this> in more convenient ways
#define self (*this)

struct MatrixSizeError : std::exception
{
private:
	const char *msg;

public:
	MatrixSizeError() : msg("Matrix Size Error") {}
	MatrixSizeError(const char *_msg) : msg(_msg) {}

	virtual const char *what() const override { return msg; }
};

template<typename T>
class Matrix
{
public: // -- enums / etc -- //



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

	Matrix(const Matrix &other) = default;
	// constructs from another matrix
	// other is guaranteed to be empty() afterwards
	Matrix(Matrix &&other) : r(other.r), c(other.c), data(std::move(other.data))
	{
		// empty other matrix
		other.r = other.c = 0;
	}

	Matrix &operator=(const Matrix &other) = default;
	// copies from another matrix via move semantics
	// equivalent to swapping the contents of the matricies
	Matrix &operator=(Matrix &&other)
	{
		using std::swap; // ADL idiom

		// mov asgn swap idiom
		swap(r, other.r);
		swap(c, other.c);
		swap(data, other.data);
	}

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
	// empties the matrix
	bool clear()
	{
		r = c = 0;
		data.clear();
	}

	// resizes the matrix to the specified dimensions
	// resizing to 0xn or nx0 is equivalent to calling clear()
	void resize(std::size_t rows, std::size_t cols)
	{
		// if either dimension was zero, clear matrix
		if (rows == 0 || cols == 0) clear();
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
			swap(self(a, i), self(b, i));
	}

	// multiplies the elements in a row by a scalar
	void multRow(std::size_t row, T f)
	{
		for (std::size_t i = 0; i < c; ++i) self(row, i) *= f;
	}
	// divides the elements in a row by a scalar
	void divRow(std::size_t row, T f)
	{
		for (std::size_t i = 0; i < c; ++i) self(row, i) /= f;
	}

	// adds row <from> to row <to>
	void addRow(std::size_t to, std::size_t from)
	{
		for (std::size_t i = 0; i < c; ++i) self(to, i) += self(from, i);
	}
	// subtracts row <from> from row <to>
	void subRow(std::size_t to, std::size_t from)
	{
		for (std::size_t i = 0; i < c; ++i) self(to, i) -= self(from, i);
	}

	// adds a multiple of row <from> to row <to>
	void addRowMult(std::size_t to, std::size_t from, T f)
	{
		for (std::size_t i = 0; i < c; ++i) self(to, i) += self(from, i) * f;
	}
	// subtracts a multiple of row <from> to row <to>
	void subRowMult(std::size_t to, std::size_t from, T f)
	{
		for (std::size_t i = 0; i < c; ++i) self(to, i) -= self(from, i) * f;
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
		case 1: return self(0, 0);
		case 2: return self(0, 0) * self(1, 1) - self(1, 0)*self(0, 1);

		default:
			T sum = 0; // initialize sum to zero

			// add cofactors along first row
			for (std::size_t col = 0; col < c; ++col)
				sum += cofactor(0, col);

			return sum;
		}
	}

	// returns the <row><col> submatrix (i.e. the result of removing row <row> and col <col>)
	// throws MatrixSizeError if matrix is empty
	Matrix submatrix(std::size_t row, std::size_t col) const
	{
		if (r == 0) throw MatrixSizeError("Cannot take a minor from an empty matrix");

		// allocate the result
		Matrix res(r - 1, c - 1);

		// fill its entries
		for (std::size_t _row = 0; _row < row; ++_row)
		{
			for (std::size_t _col = 0; _col < col; ++_col)
				res(_row, _col) = self(_row, _col);
			for (std::size_t _col = col + 1; _col < c; ++_col)
				res(_row, _col - 1) = self(_row, _col);
		}
		for (std::size_t _row = row + 1; _row < r; ++_row)
		{
			for (std::size_t _col = 0; _col < col; ++_col)
				res(_row - 1, _col) = self(_row, _col);
			for (std::size_t _col = col + 1; _col < c; ++_col)
				res(_row - 1, _col - 1) = self(_row, _col);
		}

		return res;
	}
	// returns the <row><col> minor (i.e. submatrix(row, col).det() )
	// throws by submatrix()
	T minor(std::size_t row, std::size_t col) const
	{
		return submatrix(row, col).det();
	}
	// returns the <row><col> cofactor (i.e. (-1)^(row+col) * (row,col) * minor(row,col) )
	// throws by minor()
	T cofactor(std::size_t row, std::size_t col) const
	{
		return ((row + col) & 1 ? -1 : 1) * self(row, col) * minor(row, col);
	}

	// returns the cofactor matrix of this matrix
	// throws by cofactor()
	Matrix cofactorMatrix() const
	{
		// allocate the result
		Matrix res(r, c);

		// fill with cofactors
		for (std::size_t row = 0; row < r; ++row)
			for (std::size_t col = 0; col < c; ++col)
				res(row, col) = cofactor(row, col);

		return res;
	}
	// returns the adjugate matrix (i.e. transpose of cofactor matrix)
	// throws by cofactor()
	Matrix adjugate() const
	{
		// allocate the result
		Matrix res(r, c);

		// fill with cofactors
		for (std::size_t row = 0; row < r; ++row)
			for (std::size_t col = 0; col < c; ++col)
				res(col, row) = cofactor(row, col);

		return res;
	}

	// transposes the matrix
	void transpose()
	{
		using std::swap; // ADL idiom

		// as a special case, row/column vectors will have identical flattened structures
		if (r == 1 || c == 1) swap(r, c);
		// square matricies can use simple swaps
		else if (r == c)
		{
			for (std::size_t row = 0; row < r; ++row)
				for (std::size_t col = row + 1; col < c; ++col)
					swap(self(row, col), self(col, row));
		}
		// otherwise it's comlicated. just copy a transposed version
		else *this = std::move(transposed());
	}
	// returns a copy of this matrix that has been transposed
	Matrix transposed() const
	{
		// allocate the matrix
		Matrix res(c, r);

		// copy the transposed entries
		for (std::size_t row = 0; row < r; ++row)
			for (std::size_t col = 0; col < c; ++col)
				res(col, row) = self(row, col);

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
					for (; top < bottom && self(top, col) != 0; ++top);
					// wind bottom to next nonzero entry
					for (; top < bottom && self(bottom, col) == 0; --bottom);

					// if they're still valid, perform the swap
					if (top < bottom) swapRows(top, bottom);
					// otherwise, we're done sorting
					else break;
				}

				// if after sorting this element is nonzero, it is the leading entry
				if (self(row, col) != 0) break;
			}
			// if we didn't find a leading entry, the rest of the matrix is zeroes (so we're done)
			if (col == c) return row;

			// make this row's leading entry a 1 via row division //

			// store lead entry for efficiency
			T temp = self(row, col);
			// all elements to left of leading entry are zeros, so we can ignore them
			for (std::size_t j = col + 1; j < c; ++j) self(row, j) /= temp;
			// setting leading entry to 1 is more efficient and prevents rounding errors
			self(row, col) = 1;

			// use this row to eliminate the lower rows //

			for (std::size_t j = row + 1; j < r && self(j, col) != 0; ++j)
			{
				// store factor to subtract by for efficiency
				temp = self(j, col);
				// all source elements to left of col are zero, so we can ignore them
				for (std::size_t k = col + 1; k < c; ++k) self(j, k) -= self(row, k) * temp;
				// setting entry to 0 is more efficient and prevents rounding errors
				self(j, col) = 0;
			}
		}
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
					for (; top < bottom && self(top, col) != 0; ++top);
					// wind bottom to next nonzero entry
					for (; top < bottom && self(bottom, col) == 0; --bottom);
					
					// if they're still valid, perform the swap
					if (top < bottom) swapRows(top, bottom);
					// otherwise, we're done sorting
					else break;
				}

				// if after sorting this element is nonzero, it is the leading entry
				if (self(row, col) != 0) break;
			}
			// if we didn't find a leading entry, the rest of the matrix is zeroes (so we're done)
			if (col == c) return row;

			// make this row's leading entry a 1 via row division //

			// store lead entry for efficiency
			T temp = self(row, col);
			// all elements to left of leading entry are zeros, so we can ignore them
			for (std::size_t j = col + 1; j < c; ++j) self(row, j) /= temp;
			// setting leading entry to 1 is more efficient and prevents rounding errors
			self(row, col) = 1;

			// use this row to eliminate the higher rows //

			for (std::size_t j = 0; j < row; ++j)
			{
				// store factor to subtract by for efficiency
				temp = self(j, col);
				if (temp != 0)
				{
					// all source elements to left of col are zero, so we can ignore them
					for (std::size_t k = col + 1; k < c; ++k) self(j, k) -= self(row, k) * temp;
					// setting entry to 0 is more efficient and prevents rounding errors
					self(j, col) = 0;
				}
			}

			// use this row to eliminate the lower rows //

			for (std::size_t j = row + 1; j < r && self(j, col) != 0; ++j)
			{
				// store factor to subtract by for efficiency
				temp = self(j, col);
				// all source 
				for (std::size_t k = col + 1; k < c; ++k) self(j, k) -= self(row, k) * temp;
				// setting entry to 0 is more efficient and prevents rounding errors
				self(j, col) = 0;
			}
		}

		// successfully converted to RREF
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
				util(row, col) = self(row, col);
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
Matrix<T> &operator+=(Matrix<T> &a, const Matrix<T> &b)
{
	// ensure sizes are ok
	if (a.r != b.r || a.c != b.c) throw MatrixSizeError("Matrix addition requires the matricies be of same size");

	for (std::size_t row = 0; row < r; ++row)
		for (std::size_t col = 0; col < c; ++col)
			a(row, col) += b(row, col);

	return a;
}
template<typename T>
Matrix<T> operator+(const Matrix<T> &a, const Matrix<T> &b)
{
	Matrix<T> res = a;
	return res += b;
}

// subtracts matrix b from matrix a
// throws MatrixSizeError if matricies are of different sizes
template<typename T>
Matrix<T> &operator-=(Matrix<T> &a, const Matrix<T> &b)
{
	// ensure sizes are ok
	if (a.r != b.r || a.c != b.c) throw MatrixSizeError("Matrix subtraction requires the matricies be of same size");

	for (std::size_t row = 0; row < r; ++row)
		for (std::size_t col = 0; col < c; ++col)
			a(row, col) -= b(row, col);

	return a;
}
template<typename T>
Matrix<T> operator-(const Matrix<T> &a, const Matrix<T> &b)
{
	Matrix<T> res = a;
	return res -= b;
}

// multiplies the matrix by a scalar
template<typename T>
Matrix<T> &operator*=(Matrix<T> &m, const T &f)
{
	for (std::size_t row = 0; row < r; ++row)
		for (std::size_t col = 0; col < c; ++col)
			a(row, col) *= f;
}
template<typename T>
Matrix<T> operator*(const Matrix<T> &m, const T &f)
{
	Matrix<T> res = m;
	return res *= f;
}
template<typename T>
Matrix<T> operator*(const T &f, const Matrix<T> &m)
{
	Matrix<T> res = m;
	return res *= f;
}

// multiplies matrix lhs by matrix rhs
// throws MatrixSizeError if matricies are incompatible
template<typename T>
Matrix<T> operator*(const Matrix<T> &lhs, const Matrix<T> &rhs)
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
Matrix<T> operator*=(Matrix<T> &lhs, const Matrix<T> &rhs)
{
	return *this = lhs * rhs;
}

// -- cmp operator definitions -- //

// returns true if matricies are of same size and have identical contents
template<typename T>
bool operator==(const Matrix<T> &a, const Matrix<T> &b)
{
	// different sizes are unequal by definition
	if (a.rows() != b.rows() || a.cols() != b.cols()) return false;

	// otherwise must contain identical elements
	for (std::size_t row = 0row < a.rows(); ++row)
		for (std::size_t col = 0; col < a.cols(); ++col)
			if (a(row, col) != b(row, col)) return false;

	return true;
}
template<typename T>
bool operator!=(const Matrix<T> &a, const Matrix<T> &b) { return !(a == b); }

#endif
