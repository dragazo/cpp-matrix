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

#define __MATRIX_DIAGNOSTICS 0

struct MatrixSizeError : std::exception
{
private:
	const char *msg;

public:
	MatrixSizeError() : msg("Matrix Size Error") {}
	MatrixSizeError(const char *_msg) : msg(_msg) {}

	virtual const char *what() const override { return msg; }
};

// represents a mathematical matrix
// T is assumed to be a POD value-type
// T is required to be explicitly constructable from int
// T is required to have binary operators +, -, *, / (and their compound assignments), ==, and != defined, capable of accepting T as both arguments 
// T is required to have unary operator - defined, capable of accepting T as its argument
template<typename T>
class Matrix
{
public: // -- enums / etc -- //



private: // -- data -- //

	T *data;          // the elements in the array
	std::size_t cap;  // capacity of array
	std::size_t r, c; // number of rows / cols

private: // -- helpers -- //

	// performs a generalized matrix reduction via elementary row operations, used by (for example) REF, RREF, and det
	// also returns the determinant of the matrix if square, or of the square matrix resulting from truncating a rectangular matrix to the largest square matrix it can contain
	// <kill_upper> specifies if the reduction should also eliminate the upper triangle (e.g. RREF)
	// <det>        the location to store the calculated determinant, or null to ignore
	// <only_det>   specifies that we're only interested in the determinant (i.e. will exit early if found to be zero, in which case rank may not be correct)
	// returns      the rank of the matrix (unless <only_det> triggered an early exit)
	// complexity   O(n^2)
	std::size_t reduce(bool kill_upper, T *det, bool only_det)
	{
		using std::swap; // ADL idiom

		T _det = (T)1;        // initialize the resulting determinant
		T temp;               // storage location for leading entry

		// iterate through each row
		for (std::size_t row = 0, col = 0; row < r; ++row, ++col)
		{
			// search for a leading entry
			for (; col < c && self(row, col) == (T)0; ++col)
			{
				// look down the column for a row to swap in
				for(std::size_t i = row + 1; i < r; ++i)
					if (self(i, col) != (T)0)
					{
						// swap this row in
						for (std::size_t j = col; j < c; ++j) swap(self(row, j), self(i, j));
						// negate det due to a row swap
						_det = -_det;
						// and we found our leading entry
						goto found_entry;
					}

				// otherwise we're missing a leading entry in the main diagonal, which means det is zero
				// setting col to c makes it think there was a row of zeroes, and thus det of zero
				if (only_det) col = c;
			}
			found_entry:

			// if we didn't find a leading entry, the rest of the matrix is zeroes (so we're done)
			if (col == c)
			{
				if (det) *det = (T)0; // no leading entry means a row of all zero, which means a det of zero
				return row;
			}

			// make this row's leading entry a 1 via row division //

			// store lead entry for efficiency
			temp = self(row, col);
			// all elements to left of leading entry are zeros, so we can ignore them
			for (std::size_t i = col + 1; i < c; ++i) self(row, i) /= temp;
			// setting leading entry to 1 is more efficient and prevents rounding errors
			self(row, col) = (T)1;
			_det *= temp; // account for this entry in the determinant

			// use this row to eliminate the higher rows //

			// only do this if requested
			if (kill_upper)
			{
				for (std::size_t i = 0; i < row; ++i)
				{
					// store factor to subtract by for efficiency
					temp = self(i, col);
					if (temp != (T)0)
					{
						// all source elements to left of col are zero, so we can ignore them
						for (std::size_t j = col + 1; j < c; ++j) self(i, j) -= self(row, j) * temp;
						// setting entry to 0 is more efficient and prevents rounding errors
						self(i, col) = (T)0;
					}
				}
			}

			// use this row to eliminate the lower rows //
			
			for (std::size_t i = row + 1; i < r; ++i)
			{
				// store factor to subtract by for efficiency
				temp = self(i, col);
				if (temp != (T)0)
				{
					// all source elements to left of col are zero, so we can ignore them
					for (std::size_t j = col + 1; j < c; ++j) self(i, j) -= self(row, j) * temp;
					// setting entry to 0 is more efficient and prevents rounding errors
					self(i, col) = (T)0;
				}
			}
		}

		// export the calculated determinant
		if (det) *det = _det;

		// successfully converted to RREF
		return r;
	}

public: // -- ctor / dtor / asgn -- //

	// creates an empty matrix
	Matrix() : data(nullptr), cap(0), r(0), c(0)
	{
		#if __MATRIX_DIAGNOSTICS
		std::cout << "def ctor\n";
		#endif
	}

	// creates a matrix with the specified dimensions
	// element contents are undefined
	Matrix(std::size_t rows, std::size_t cols)
	{
		// if either dimension was zero, both are zero
		if (rows == 0 || cols == 0)
		{
			cap = r = c = 0;
			data = nullptr;
		}
		else
		{
			cap = rows * cols;
			data = new T[cap];
			r = rows;
			c = cols;
		}

		#if __MATRIX_DIAGNOSTICS
		std::cout << "arg ctor\n";
		#endif
	}

	~Matrix()
	{
		// free the array
		delete[] data;
	}

	Matrix(const Matrix &other) : r(other.r), c(other.c)
	{
		// allocate and copy the bare minimum
		cap = r * c;
		data = cap != 0 ? new T[cap] : nullptr;
		for (std::size_t i = 0; i < cap; ++i) data[i] = other.data[i];

		#if __MATRIX_DIAGNOSTICS
		std::cout << "cpy ctor\n";
		#endif
	}
	// constructs from another matrix by using its allocated resources
	// other is guaranteed to be empty() afterwards
	Matrix(Matrix &&other) : data(other.data), cap(other.cap), r(other.r), c(other.c)
	{
		// empty other matrix
		other.cap = other.r = other.c = 0;
		other.data = nullptr;

		#if __MATRIX_DIAGNOSTICS
		std::cout << "mov ctor\n";
		#endif
	}

	Matrix &operator=(const Matrix &other)
	{
		// size up to fit other
		resize_dump(other.r, other.c);
		// copy the data over
		for (std::size_t i = 0; i < r * c; ++i) data[i] = other.data[i];

		#if __MATRIX_DIAGNOSTICS
		std::cout << "cpy asgn\n";
		#endif

		return *this;
	}
	// copies from another matrix via move semantics
	// equivalent to swapping the contents of the matricies
	Matrix &operator=(Matrix &&other)
	{
		using std::swap; // ADL idiom

		// mov asgn swap idiom
		swap(data, other.data);
		swap(cap, other.cap);
		swap(r, other.r);
		swap(c, other.c);

		#if __MATRIX_DIAGNOSTICS
		std::cout << "mov asgn\n";
		#endif

		return *this;
	}

public: // -- utilities -- //

	// gets the number of rows/cols in the matrix
	std::size_t rows() const { return r; }
	std::size_t cols() const { return c; }

	// returns the number of elements in the matrix
	std::size_t size() const { return r * c; }
	// returns the current capacity of the matrix
	std::size_t capacity() const { return cap; }

	// returns true if the matrix is square (including empty)
	bool square() const { return r == c; }

	// returns true if the matrix is empty
	bool empty() const { return r == 0; }
	// sets the matrix to the "empty" state
	void clear() { r = c = 0; }

	// resizes the matrix to the specified dimensions, not making any attempt to preserve the contents
	// the contents of the result are undefined except that resizing to nx0 or 0xn is equivalent to clear() and resizing to same size is no-op
	void resize_dump(std::size_t rows, std::size_t cols)
	{
		r = rows;
		c = cols;

		// reallocate if we don't have enough space
		if (rows * cols > cap)
		{
			cap = rows * cols;
			delete[] data; // delete on null is defined to be no-op
			data = new T[cap];
		}
	}
	// resizes the matrix to the specified dimensions, preserving the contents after the call
	// this ensures the i,j elements before and after are equal over the region in which both sizes were defined
	// reducing a dimension truncates those values and expanded sections are undefined. no change is no-op
	void resize(std::size_t rows, std::size_t cols)
	{
		// save previous size
		std::size_t _rows = r, _cols = c;
		// get the smallest size values
		std::size_t min_r = (_rows < rows ? _rows : rows), min_c = (_cols < cols ? _cols : cols);

		// apply the size change
		r = rows;
		c = cols;
		reserve(rows * cols);

		// because we're using a flattened array, we only need to do some sewing if we changed #cols and it's not a row vector conversion
		// also weed out the degenerate case of resizing to/from empty
		if (cols != _cols && (_rows != 1 || rows != 1) && (min_r != 0 && min_c != 0))
		{
			// fix up the data
			for (std::size_t from = (min_r - 1) * _cols, to = (min_r - 1) * cols; from >= _cols; from -= _cols, to -= cols)
			{
				// copy row entries starting at <from> to their new locations starting at <to>
				for (std::size_t i = 0; i < min_c; ++i) data[to + i] = data[from + i];
			}
		}
	}

	// requests the matrix to set aside space for the specified number of elements
	// contents of the matrix are preserved
	void reserve(std::size_t count)
	{
		// if we don't have enough space
		if (count > cap)
		{
			// generate the new array
			cap = count;
			T *newdata = new T[cap];
			for (std::size_t i = 0; i < r * c; ++i) newdata[i] = data[i];

			// replace old array
			delete[] data;
			data = newdata;
		}
	}
	// requests the matrix to remove unused allocated space
	// contents of the matrix are preserved
	void shrink_to_fit()
	{
		// if we have too much space
		if (r * c < cap)
		{
			// generate the new array
			cap = r * c;
			T *newdata = new T[cap];
			for (std::size_t i = 0; i < cap; ++i) newdata[i] = data[i];

			// replace old array
			delete[] data;
			data = newdata;
		}
	}

	// returns the item at the specified index
	T &operator[](std::size_t index) { return data[index]; }
	T operator[](std::size_t index) const { return data[index]; }

	// gets the element at the specified index, but with bounds checking
	// throws std::out_of_range if element is out of bounds
	T &ati(std::size_t index)
	{
		if (index >= r * c) throw std::out_of_range("Specified matrix element out of bounds");
		return data[index];
	}
	T ati(std::size_t index) const
	{
		if (index >= r * c) throw std::out_of_range("Specified matrix element out of bounds");
		return data[index];
	}

	// returns the index of the specified row and col as if in a flattened array
	std::size_t index(std::size_t row, std::size_t col) const { return row * c + col; }
	// gets the row and col of the specified index
	void rowcol(std::size_t index, std::size_t &row, std::size_t &col) { row = index / c; col = index % c; }

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

	// concatenates <other> onto the bottom of the matrix
	// throw MatrixSizeError if different #cols
	void cat_rows(const Matrix &other)
	{
		if (c != other.c) throw MatrixSizeError("Cannot vertically concatenate matricies with different #cols");

		std::size_t _r = r;     // store previous row count (also puts it on the stack, so faster access)
		resize(r + other.r, c); // resize to accommodate other

		// add data from other
		for (std::size_t row = 0; row < other.r; ++row)
			for (std::size_t col = 0; col < other.c; ++col)
				self(_r + row, col) = other(row, col);
	}
	// concatenates <other> onto the right of the matrix
	// throws MatrixSizeError if matricies have different #rows
	void cat_cols(const Matrix &other)
	{
		if (r != other.r) throw MatrixSizeError("Cannot horizontally concatenate matricies with different #rows");

		std::size_t _c = c;     // store previous col count (also puts it on the stack, so faster access)
		resize(r, c + other.c); // resize to accommodate other

		// add data from other
		for (std::size_t row = 0; row < other.r; ++row)
			for (std::size_t col = 0; col < other.c; ++col)
				self(row, _c + col) = other(row, col);
	}

	// concatenates the matricies vertically
	// throws by cat_rows() member func
	static Matrix cat_rows(const Matrix &top, const Matrix &bottom) { Matrix res = top; res.cat_rows(bottom); return res; }
	static Matrix &&cat_rows(Matrix &&top, const Matrix &bottom) { top.cat_rows(bottom); return std::move(top); }
	// concatenates the matricies horizontally
	// throws by cat_cols() member func
	static Matrix cat_cols(const Matrix &left, const Matrix &right) { Matrix res = left; left.cat_cols(right); return res; }
	static Matrix &&cat_cols(Matrix &&left, const Matrix &right) { left.cat_cols(right); return std::move(top); }

public: // -- elementary row operations -- //

	// swaps rows a and b
	void swapRows(std::size_t a, std::size_t b)
	{
		using std::swap; // ADL idiom

		// early exit if swapping to same destination (one if statement potentially prevents c unnecessary swaps)
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
	// complexity: O(n^2)
	T det() const &
	{
		if (r != c) throw MatrixSizeError("Only square matricies have a determinant");
		if (r == 0) throw MatrixSizeError("Cannot take the determinant of an empty matrix");
		
		// use our handy dandy reduce function on a copy to find the determinant fancily
		Matrix cpy = self;
		T res;
		cpy.reduce(false, &res, true);
		return res;
	}
	T det() &&
	{
		if (r != c) throw MatrixSizeError("Only square matricies have a determinant");
		if (r == 0) throw MatrixSizeError("Cannot take the determinant of an empty matrix");

		// use our handy dandy reduce function to find the determinant fancily
		T res;
		reduce(false, &res, true);
		return res;
	}

	// returns true if the matrix is invertible
	bool invertible() const & { return r != 0 && r == c && det() != (T)0; }
	bool invertible() && { return r != 0 && r == c && std::move(self).det() != (T)0; }

	// returns the <row><col> submatrix (i.e. the result of removing row <row> and col <col>)
	// throws MatrixSizeError if matrix is empty
	Matrix submatrix(std::size_t row, std::size_t col) const
	{
		if (r == 0) throw MatrixSizeError("Cannot take a submatrix from an empty matrix");

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
	// throws by submatrix() and det()
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
	// throws MatrixSizeError if non-square matrix or empty
	// also throws by cofactor()
	// complexity: O(n^2) if invertible, otherwise O(n^4)
	Matrix cofactorMatrix() const
	{
		if (r != c) throw MatrixSizeError("Cannot take cofactor of non-square matrix");
		if (r == 0) throw MatrixSizeError("Cannot take cofactor of empty matrix");

		// allocate space for result
		Matrix res(r, r);
		T _det;

		// definition of adjugate: A * adj(A) = det(A) * I
		// if A is invertible:     adj(A) = det(A) * inv(A)
		// then, since the adjugate is transpose of the cofactor matrix, just transpose adjugate
		if (inverse(res, &_det))
		{
			res *= _det;
			res.transpose();
		}
		// otherwise we must resort to the other definition: adj(A) = cofactorMatrix(A).transpose()
		// then, since the adjugate is transpose of the cofactor matrix, just transpose adjugate
		else
		{
			// fill with cofactors
			for (std::size_t row = 0; row < r; ++row)
				for (std::size_t col = 0; col < c; ++col)
					res(row, col) = cofactor(row, col);
		}

		return res;
	}
	// returns the adjugate matrix (i.e. transpose of cofactor matrix)
	// throws MatrixSizeError if non-square matrix or empty
	// also throws by cofactor()
	// complexity: O(n^2) if invertible, otherwise O(n^4)
	Matrix adjugate() const
	{
		if (r != c) throw MatrixSizeError("Cannot take adjugate of non-square matrix");
		if (r == 0) throw MatrixSizeError("Cannot take adjugate of empty matrix");

		// allocate space for result
		Matrix res(r, r);
		T _det;

		// definition of adjugate: A * adj(A) = det(A) * I
		// if A is invertible:     adj(A) = det(A) * inv(A)
		if (inverse(res, &_det)) res *= _det;
		// otherwise we must resort to the other definition: adj(A) = cofactorMatrix(A).transpose()
		else
		{
			// fill with cofactors
			for (std::size_t row = 0; row < r; ++row)
				for (std::size_t col = 0; col < c; ++col)
					res(col, row) = cofactor(row, col);
		}

		return res;
	}

	// CHECK FOR OPTIMIZATIONS IN TRANSPOSE:

	// transposes the matrix
	Matrix &transpose()
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
		// otherwise it's complicated. just copy a transposed version
		else self = getTranspose();

		return self;
	}
	// returns a copy of this matrix that has been transposed
	Matrix getTranspose() const
	{
		// allocate the matrix
		Matrix res(c, r);

		// copy the transposed entries
		for (std::size_t row = 0; row < r; ++row)
			for (std::size_t col = 0; col < c; ++col)
				res(col, row) = self(row, col);

		return res;
	}

	// extracts row <row> from the matrix (as a row vector)
	Matrix getrow(std::size_t row) const
	{
		Matrix res(1, c);

		for (std::size_t i = 0; i < c; ++i) res(0, i) = self(row, i);

		return res;
	}
	// extracts column <col> from the matrix (as a column vector)
	Matrix getcol(std::size_t col) const
	{
		Matrix res(r, 1);

		for (std::size_t i = 0; i < r; ++i) res(i, 0) = self(i, col);

		return res;
	}

	// puts the matrix into row echelon form. returns the rank of the matrix
	// optionally also returns the determinant of the largest square matrix that this matrix can contain
	// complexity: O(n^2)
	std::size_t REF(T *det = nullptr) { return reduce(false, det, false); }
	// puts the matrix into reduced row echelon form. returns the rank of the matrix
	// complexity: O(n^2)
	std::size_t RREF(T *det = nullptr) { return reduce(true, det, false); }

	// attempts to find the inverse of the matrix. if successful, stores inverse in <dest> and returns true
	// optionally also returns the determinant of the largest square matrix that this matrix can contain
	// throws MatrixSizeError if matrix is not square or is empty
	// complexity: O(n^2)
	bool inverse(Matrix &dest, T *det = nullptr) const
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
		// short-circuit on zero det, as this means it's not invertible anyway
		if (util.reduce(true, det, true) == r)
		{
			// allocate dest
			dest.resize_dump(r, c);

			// copy inverse (right side) to dest
			for (std::size_t row = 0; row < r; ++row)
				for (std::size_t col = 0; col < c; ++col)
					dest(row, col) = util(row, col + c);

			return true;
		}
		else return false;
	}
	// attempts to invert the matrix. if successful, stores inverse in this matrix and returns true
	// optionally also returns the determinant of the largest square matrix that this matrix can contain
	// equivalent to m.inverse(m, det)
	// throws by inverse()
	bool invert(T *det = nullptr) { return inverse(self, det); }
};

// -- io -- //

// outputs the matrix to the stream
template<typename T> std::ostream &operator<<(std::ostream &ostr, const Matrix<T> &m)
{
	for (std::size_t row = 0; row < m.rows(); ++row)
	{
		for (std::size_t col = 0; col < m.cols(); ++col)
			ostr << std::setw(10) << m(row, col) << ' ';

		ostr << '\n';
	}

	return ostr;
}
// reads items to fill a matrix (uses the matrix's current size)
template<typename T> std::istream &operator>>(std::istream &istr, Matrix<T> &m)
{
	for (std::size_t row = 0; row < m.rows(); ++row)
		for (std::size_t col = 0; col < m.cols(); ++col)
			istr >> m(row, col);

	return istr;
}

// -- operator definitions -- //

// adds matrix b to matrix a
// throws MatrixSizeError if matricies are of different sizes
template<typename T> Matrix<T> &operator+=(Matrix<T> &a, const Matrix<T> &b)
{
	// ensure sizes are ok
	if (a.rows() != b.rows() || a.cols() != b.cols()) throw MatrixSizeError("Matrix addition requires the matricies be of same size");

	for (std::size_t row = 0; row < a.rows(); ++row)
		for (std::size_t col = 0; col < a.cols(); ++col)
			a(row, col) += b(row, col);

	return a;
}
template<typename T> Matrix<T> operator+(const Matrix<T> &a, const Matrix<T> &b) { Matrix<T> res = a; res += b; return res; }
template<typename T> Matrix<T> &&operator+(Matrix<T> &&a, const Matrix<T> &b) { a += b; return std::move(a); }
template<typename T> Matrix<T> &&operator+(const Matrix<T> &a, Matrix<T> &&b) { b += a; return std::move(b); }
template<typename T> Matrix<T> &&operator+(Matrix<T> &&a, Matrix<T> &&b) { a += b; return std::move(a); }

// subtracts matrix b from matrix a
// throws MatrixSizeError if matricies are of different sizes
template<typename T> Matrix<T> &operator-=(Matrix<T> &a, const Matrix<T> &b)
{
	// ensure sizes are ok
	if (a.rows() != b.rows() || a.cols() != b.cols()) throw MatrixSizeError("Matrix subtraction requires the matricies be of same size");

	for (std::size_t row = 0; row < a.rows(); ++row)
		for (std::size_t col = 0; col < a.cols(); ++col)
			a(row, col) -= b(row, col);

	return a;
}
template<typename T> Matrix<T> operator-(const Matrix<T> &a, const Matrix<T> &b) { Matrix<T> res = a; res -= b; return res; }
template<typename T> Matrix<T> &&operator-(Matrix<T> &&a, const Matrix<T> &b) { a -= b; return std::move(a); }
template<typename T> Matrix<T> &&operator-(const Matrix<T> &a, Matrix<T> &&b) { b -= a; return std::move(b); }
template<typename T> Matrix<T> &&operator-(Matrix<T> &&a, Matrix<T> &&b) { a -= b; return std::move(a); }

// multiplies the matrix by a scalar
template<typename T> Matrix<T> &operator*=(Matrix<T> &m, T f)
{
	for (std::size_t row = 0; row < m.rows(); ++row)
		for (std::size_t col = 0; col < m.cols(); ++col)
			m(row, col) *= f;

	return m;
}
template<typename T> Matrix<T> operator*(const Matrix<T> &m, T f) { Matrix<T> res = m; res *= f; return res; }
template<typename T> Matrix<T> operator*(T f, const Matrix<T> &m) { Matrix<T> res = m; res *= f; return res; }
template<typename T> Matrix<T> &&operator*(Matrix<T> &&m, T f) { m *= f; return std::move(m); }
template<typename T> Matrix<T> &&operator*(T f, Matrix<T> &&m) { m *= f; return std::move(m); }

// multiplies matrix lhs by matrix rhs
// throws MatrixSizeError if matricies are incompatible
template<typename T> Matrix<T> operator*(const Matrix<T> &lhs, const Matrix<T> &rhs)
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
template<typename T> Matrix<T> &operator*=(Matrix<T> &lhs, const Matrix<T> &rhs) { return lhs = lhs * rhs; }

// -- cmp operator definitions -- //

// returns true if matricies are of same size and have identical contents
template<typename T> bool operator==(const Matrix<T> &a, const Matrix<T> &b)
{
	// different sizes are unequal by definition
	if (a.rows() != b.rows() || a.cols() != b.cols()) return false;

	// otherwise must contain identical elements
	for (std::size_t row = 0; row < a.rows(); ++row)
		for (std::size_t col = 0; col < a.cols(); ++col)
			if (a(row, col) != b(row, col)) return false;

	return true;
}
template<typename T> bool operator!=(const Matrix<T> &a, const Matrix<T> &b) { return !(a == b); }

#endif
