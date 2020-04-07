#ifndef DRAGAZO_MATRICES_H
#define DRAGAZO_MATRICES_H

#include <utility>
#include <type_traits>
#include <vector>
#include <exception>
#include <algorithm>
#include <execution>
#include <iterator>
#include <functional>

namespace matrix_impl
{
	struct matrix_size_error : std::exception { using std::exception::exception; };

	template<typename T, typename D = T>
	struct val_iter
	{
		typedef std::random_access_iterator_tag iterator_category;

		typedef D    difference_type;
		typedef void value_type;
		typedef T    reference;
		typedef void pointer;

		static_assert(std::is_integral_v<T>, "T must be an integral type");
		static_assert(std::is_integral_v<D>, "D must be an integral type");

		T val;

		T operator*() const noexcept { return val; }

		val_iter &operator++() noexcept { ++val; return *this; }
		val_iter operator++(int) noexcept { return { val++ }; }

		val_iter &operator--() noexcept { --val; return *this; }
		val_iter operator--(int) noexcept { return { val-- }; }

		val_iter &operator+=(D off) noexcept { val += static_cast<T>(off); return *this; }
		friend val_iter operator+(val_iter i, D off) noexcept { i += off; return i; }
		friend val_iter operator+(D off, val_iter i) noexcept { i += off; return i; }

		val_iter &operator-=(D off) noexcept { val -= static_cast<T>(off); return *this; }
		friend val_iter operator-(val_iter i, D off) noexcept { i -= off; return i; }

		friend D operator-(val_iter a, val_iter b) noexcept { return static_cast<D>(a.val - b.val); }

		friend bool operator==(val_iter a, val_iter b) noexcept { return a.val == b.val; }
		friend bool operator!=(val_iter a, val_iter b) noexcept { return a.val != b.val; }
	};

	struct noop
	{
		template<typename A>
		decltype(auto) operator()(A &&a) const noexcept { return std::forward<A>(a); }
	};
	struct mover
	{
		template<typename A>
		decltype(auto) operator()(A &a) const noexcept { return std::move(a); }
	};

	// ----------------------------------------------------------------------------------------------------

	template<typename A, typename OP>
	struct unary_row_expr
	{
	private:
		const A &a;
	public:
		static_assert(std::is_same_v<std::remove_cv_t<A>, A>, "A must not be cv qualified");
		unary_row_expr(const A &_a) noexcept : a(_a) {}
		decltype(auto) operator[](std::size_t i) const { return OP{}(a[i]); }
		std::size_t size() const { return a.size(); }
	};
	template<typename A, typename B, typename OP>
	struct binary_row_expr
	{
	private:
		const A &a;
		const B &b;
	public:
		static_assert(std::is_same_v<std::remove_cv_t<A>, A>, "A must not be cv qualified");
		static_assert(std::is_same_v<std::remove_cv_t<B>, B>, "B must not be cv qualified");
		binary_row_expr(const A &_a, const B &_b) : a(_a), b(_b) { assert(a.size() == b.size()); }
		decltype(auto) operator[](std::size_t i) const { return OP{}(a[i], b[i]); }
		std::size_t size() const { return a.size(); }
	};
	template<typename Vec, typename Scalar, typename OP, bool left_scalar>
	struct row_scalar_expr
	{
	private:
		const Vec &vec;
		Scalar     scalar;
	public:
		static_assert(std::is_same_v<std::remove_cv_t<Vec>, Vec>, "Vec must not be cv qualified");
		static_assert(std::is_same_v<std::remove_cv_t<std::remove_reference_t<Scalar>>, Scalar>, "Scalar must not be ref or cv qualified");
		row_scalar_expr(const Vec &v, const Scalar &s) : vec(v), scalar(s) {}
		row_scalar_expr(const Vec &v, Scalar &&s) : vec(v), scalar(std::move(s)) {}
		decltype(auto) operator[](std::size_t i) const
		{
			if constexpr (left_scalar) return OP{}(scalar, vec[i]);
			else return OP{}(vec[i], scalar);
		}
		std::size_t size() const { return vec.size(); }
	};

	// ----------------------------------------------------------------------------------------------------

	template<typename T> struct is_row_expr : std::false_type {};
	template<typename T> struct is_row_view : std::false_type {};

	template<typename T> inline constexpr bool experizable = is_row_expr<T>::value || is_row_view<T>::value;

	template<typename T> struct is_unseq : std::bool_constant<std::is_integral_v<T> || std::is_floating_point_v<T>> {};
	template<typename T> struct is_unseq<T&> : std::bool_constant<is_unseq<T>::value> {};
	template<typename T> struct is_unseq<T&&> : std::bool_constant<is_unseq<T>::value> {};
	template<typename T> struct is_unseq<const T> : std::bool_constant<is_unseq<T>::value> {};

	template<typename ...T>
	auto &exec_policy()
	{
		static_assert(sizeof...(T) > 0);
		if constexpr ((... && is_unseq<T>::value)) return std::execution::par_unseq;
		else return std::execution::par;
	}

	template<typename T, typename Allocator> class matrix;

	template<typename T_cv> struct _elem_iter;
	template<typename T_cv> struct _row_view;
	template<typename T_cv> struct _row_iter;

	// element iterator suffices to be a raw pointer, but wrap it to make them distinct types
	template<typename T_cv>
	struct _elem_iter
	{
	private: // -- data -- //

		T_cv *p;

		friend struct _row_view<T_cv>;

		_elem_iter(T_cv *_p) noexcept : p(_p) {}

	public: // -- interface -- //

		_elem_iter() noexcept = default;

		typedef std::random_access_iterator_tag iterator_category;

		typedef std::ptrdiff_t difference_type;
		typedef T_cv           value_type;
		typedef T_cv          &reference;
		typedef T_cv          *pointer;

		reference operator*() const noexcept { return *p; }

		_elem_iter &operator++() noexcept { ++p; return *this; }
		_elem_iter operator++(int) noexcept { return { p++ }; }

		_elem_iter &operator--() noexcept { --p; return *this; }
		_elem_iter operator--(int) noexcept { return { p-- }; }

		_elem_iter &operator+=(difference_type d) noexcept { p += d; return *this; }
		friend _elem_iter operator+(_elem_iter i, difference_type d) noexcept { return { i.p + d }; }
		friend _elem_iter operator+(difference_type d, _elem_iter i) noexcept { return { i.p + d }; }

		_elem_iter &operator-=(difference_type d) noexcept { p -= d; return *this; }
		friend _elem_iter operator-(_elem_iter i, difference_type d) noexcept { return { i.p - d }; }

		friend difference_type operator-(_elem_iter a, _elem_iter b) noexcept { return a.p - b.p; }

		friend bool operator==(_elem_iter a, _elem_iter b) noexcept { return a.p == b.p; }
		friend bool operator!=(_elem_iter a, _elem_iter b) noexcept { return a.p != b.p; }
	};

	// represents a non-owning view into a row. assignment to this object just repoints the view.
	template<typename T_cv>
	struct _row_view
	{
	private: // -- data -- //

		T_cv       *vals;
		std::size_t count;

		template<typename A, typename B> friend class matrix;
		template<typename A> friend struct _row_iter;

		_row_view(T_cv *v, std::size_t c) noexcept : vals(v), count(c) {}

	public: // -- traits -- //

		typedef T_cv  value_type;
		typedef T_cv &reference;

		typedef std::size_t    size_type;
		typedef std::ptrdiff_t difference_type;

		typedef _elem_iter<T_cv>                iterator;
		typedef std::reverse_iterator<iterator> reverse_iterator;

	public: // -- interface -- //

		// constructs a view object in an unusable state (must be properly initialized before use)
		_row_view() noexcept = default;

		_row_view(const _row_view&) noexcept = default;

		// repoints this object to view a different row.
		// DOES NOT COPY THE ROW.
		_row_view &operator=(const _row_view&) & noexcept = default;
	
		// given a vector expression, computes and assigns its values to this row in parallel
		template<typename E, std::enable_if_t<is_row_expr<E>::value, int> = 0>
		const _row_view &operator=(const E &expr) const
		{
			assert(count == expr.size());
			std::for_each_n(exec_policy<T_cv&, decltype(expr[0])>(), val_iter<std::size_t>{ 0 }, count, [this, expr](std::size_t i) { vals[i] = expr[i]; });
			return *this;
		}
		template<typename E, std::enable_if_t<experizable<E>, int> = 0>
		const _row_view &operator+=(const E &expr) const
		{
			assert(count == expr.size());
			std::for_each_n(exec_policy<T_cv&, decltype(expr[0])>(), val_iter<std::size_t>{ 0 }, count, [this, expr](std::size_t i) { vals[i] += expr[i]; });
			return *this;
		}
		template<typename E, std::enable_if_t<experizable<E>, int> = 0>
		const _row_view &operator-=(const E &expr) const
		{
			assert(count == expr.size());
			std::for_each_n(exec_policy<T_cv&, decltype(expr[0])>(), val_iter<std::size_t>{ 0 }, count, [this, expr](std::size_t i) { vals[i] -= expr[i]; });
			return *this;
		}

		// for use in row expressions - effectively an xvalue view into this row
		template<int _ = 0, std::enable_if_t<_ == 0 && std::is_same_v<T_cv, std::remove_cv_t<T_cv>>, int> = 0>
		decltype(auto) move() const { return unary_row_expr<_row_view, mover>{ *this }; }

		operator _row_view<const T_cv>() const noexcept { return { vals, count }; }

		// iterates through each item in the row
		iterator begin() const noexcept { return { vals }; }
		iterator end() const noexcept { return { vals + count }; }

		// iterates through each item in the row in reverse order
		reverse_iterator rbegin() const noexcept { return { end() }; }
		reverse_iterator rend() const noexcept { return { begin() }; }

		// gets the size of this row (number of columns)
		size_type size() const noexcept { return count; }
		size_type cols() const noexcept { return count; }

		// accesses the element in the specified column (of this row)
		reference operator[](std::size_t i) const { return vals[i]; }
		reference at(std::size_t i) const { return i < count ? vals[i] : throw std::out_of_range("index out of bounds"); }

		// returns the first/last item in the row (und if row is invalid)
		reference front() const { return vals[0]; }
		reference back() const { return vals[count - 1]; }
	};

	template<typename T_cv>
	struct _row_iter
	{
	private: // -- data -- //

		T_cv       *vals;
		std::size_t count;

		template<typename T, typename Allocator> friend class matrix;

		_row_iter(T_cv *v, std::size_t c) noexcept : vals(v), count(c) {}

	public: // -- traits -- //

		typedef std::random_access_iterator_tag iterator_category;

		typedef std::ptrdiff_t  difference_type;
		typedef void            value_type;
		typedef _row_view<T_cv> reference;
		typedef void            pointer;

	public: // -- interface -- //

		// constructs an iterator in an unusable state (must be properly initialized before use)
		_row_iter() noexcept = default;

		operator _row_iter<const T_cv>() const noexcept { return { vals, count }; }

		reference operator*() const { return { vals, count }; }
		reference operator[](difference_type off) const { return { vals + off * count, count }; }

		_row_iter &operator++() noexcept { vals += count; return *this; }
		_row_iter operator++(int) noexcept { auto cpy = *this; ++*this; return cpy; }

		_row_iter &operator--() noexcept { vals -= count; return *this; }
		_row_iter operator--(int) noexcept { auto cpy = *this; --*this; return cpy; }

		_row_iter &operator+=(difference_type off) noexcept { vals += static_cast<std::ptrdiff_t>(count) * off; return *this; }
		friend _row_iter operator+(_row_iter a, difference_type off) noexcept { a += off; return a; }

		_row_iter &operator-=(difference_type off) noexcept { vals -= static_cast<std::ptrdiff_t>(count) * off; return *this; }
		friend _row_iter operator-(_row_iter a, difference_type off) noexcept { a -= off; return a; }

		friend difference_type operator-(_row_iter a, _row_iter b) noexcept { return (a.vals - b.vals) / static_cast<std::ptrdiff_t>(a.count); }

		friend bool operator==(_row_iter a, _row_iter b) noexcept { return a.vals == b.vals; }
		friend bool operator!=(_row_iter a, _row_iter b) noexcept { return a.vals != b.vals; }
	};

	// ----------------------------------------------------------------------------------------------------

	template<typename A, typename OP> struct is_row_expr<unary_row_expr<A, OP>> : std::true_type {};
	template<typename A, typename B, typename OP> struct is_row_expr<binary_row_expr<A, B, OP>> : std::true_type {};
	template<typename Vec, typename Scalar, typename OP, bool left_scalar> struct is_row_expr<row_scalar_expr<Vec, Scalar, OP, left_scalar>> : std::true_type {};

	template<typename T_cv> struct is_unseq<_row_view<T_cv>> : std::bool_constant<is_unseq<std::remove_cv_t<T_cv>>::value> {};

	template<typename A, typename OP> struct is_unseq<unary_row_expr<A, OP>> : std::bool_constant<is_unseq<A>::value> {};
	template<typename A, typename B, typename OP> struct is_unseq<binary_row_expr<A, B, OP>> : std::bool_constant<is_unseq<A>::value && is_unseq<B>::value> {};
	template<typename Vec, typename Scalar, typename OP, bool left_scalar> struct is_unseq<row_scalar_expr<Vec, Scalar, OP, left_scalar>> : std::bool_constant<is_unseq<decltype(std::declval<Vec>()[0])>::value && is_unseq<Scalar>::value> {};

	template<typename T_cv> struct is_row_view<_row_view<T_cv>> : std::true_type {};

	// ----------------------------------------------------------------------------------------------------

	template<typename A, std::enable_if_t<experizable<A>, int> = 0>
	auto operator+(const A &a) noexcept { return unary_row_expr<A, noop>{ a }; }
	template<typename A, std::enable_if_t<experizable<A>, int> = 0>
	auto operator-(const A &a) noexcept { return unary_row_expr<A, std::negate<>>{ a }; }

	template<typename A, typename B, std::enable_if_t<experizable<A> && experizable<B>, int> = 0>
	auto operator+(const A &a, const B &b) noexcept { return binary_row_expr<A, B, std::plus<>>{ a, b }; }
	template<typename A, typename B, std::enable_if_t<experizable<A> && experizable<B>, int> = 0>
	auto operator-(const A &a, const B &b) noexcept { return binary_row_expr<A, B, std::minus<>>{ a, b }; }

	template<typename Vec, typename Scalar, std::enable_if_t<experizable<Vec> && !std::is_same_v<decltype((*(const Vec*)0)[0] * (*(const std::decay_t<Scalar>*)0)), void> , int> = 0>
	auto operator*(const Vec &vec, Scalar &&scalar) noexcept { return row_scalar_expr<Vec, std::decay_t<Scalar>, std::multiplies<>, false>{ vec, std::forward<Scalar>(scalar) }; }
	template<typename Vec, typename Scalar, std::enable_if_t<experizable<Vec> && !std::is_same_v<decltype((*(const std::decay_t<Scalar>*)0) * (*(const Vec*)0)[0]), void>, int> = 0>
	auto operator*(Scalar &&scalar, const Vec &vec) noexcept { return row_scalar_expr<Vec, std::decay_t<Scalar>, std::multiplies<>, true>{ vec, std::forward<Scalar>(scalar) }; }

	template<typename Vec, typename Scalar, std::enable_if_t<experizable<Vec> && !std::is_same_v<decltype((*(const Vec*)0)[0] * (*(const std::decay_t<Scalar>*)0)), void>, int> = 0>
	auto operator/(const Vec &vec, Scalar &&scalar) noexcept { return row_scalar_expr<Vec, std::decay_t<Scalar>, std::divides<>, false>{ vec, std::forward<Scalar>(scalar) }; }

	// ----------------------------------------------------------------------------------------------------

	// represents the mathematical matrix construct with elements of type T.
	// T is required to have operators (+, -, *, /) defined (TxT -> T), as well as their compound assignments, and behave similarly to the real numbers (int is allowed but may produce incorrect results for e.g. rref or det due to truncation on divide).
	// T is required to be default constructible.
	// T is required to be constructible via the expression (T)v where v is int.
	// the above operations must be safe in parallel context.
	template<typename T, typename Allocator> //typename Allocator = vectorizable_allocator<T>>
	class matrix
	{
	public: // -- requirements -- //

		static_assert(std::is_same_v<T, std::remove_cv_t<std::remove_reference_t<T>>>, "T cannot be ref or cv qualified");
		static_assert(std::is_default_constructible_v<T>, "T is required to be default constructible");
		static_assert(std::is_constructible_v<T, int>, "T is required to be constructible from int");

	private: // -- data -- //

		std::vector<T, Allocator> data;  // the flattened array of matrix elements
		std::size_t               r = 0; // number of rows
		std::size_t               c = 0; // number of cols (if either of these is zero they must both be zero)

	public: // -- related types -- //

		typedef T         value_type;
		typedef Allocator allocator_type;

		typedef _row_view<T>       row_view;
		typedef _row_view<const T> const_row_view;

		typedef _row_view<T>       flat_view;
		typedef _row_view<const T> const_flat_view;

		typedef _row_iter<T>       row_iter;
		typedef _row_iter<const T> const_row_iter;

		typedef std::reverse_iterator<row_iter>       reverse_row_iter;
		typedef std::reverse_iterator<const_row_iter> const_reverse_row_iter;

		typedef row_iter       iterator;
		typedef const_row_iter const_iterator;

		typedef reverse_row_iter       reverse_iterator;
		typedef const_reverse_row_iter const_reverse_iter;

	private: // -- helper functions -- //

		// effectively returns |a - b| but does not require abs() to be overloaded
		static T dist(const T &a, const T &b)
		{
			return a >= b ? a - b : b - a;
		}

		// performs a generalized matrix reduction via elementary row operations, used by (for example) REF, RREF, and det.
		// also returns the determinant of the matrix if square, or of the square matrix resulting from truncating a rectangular matrix to the largest square matrix it can contain.
		// kill_upper - specifies if the reduction should also eliminate the upper triangle (e.g. RREF)
		// det        - the location to store the calculated determinant, or null to ignore
		// only_det   - specifies that we're only interested in the determinant (i.e. will exit early if found to be zero, in which case rank may not be correct)
		// returns the rank of the matrix (unless only_det triggered an early exit)
		std::size_t reduce(bool kill_upper, T *det, bool only_det)
		{
			using std::swap;

			const T zero{ 0 };

			T _det{ 1 }; // initialize the resulting determinant
			T temp;      // storage location for leading entry

			// iterate through the matrix along a diagonal
			for (std::size_t row = 0, col = 0; row < r; ++row, ++col)
			{
				// search for a leading entry
				for (; col < c && operator()(row, col) == zero; ++col)
				{
					// look down the column for a row to swap in
					for (std::size_t i = row + 1; i < r; ++i)
						if (operator()(i, col) != zero)
						{
							// swap this row in
							for (std::size_t j = col; j < c; ++j) swap(operator()(row, j), operator()(i, j));
							_det = -std::move(_det); // negate det due to a row swap
							goto found_entry; // stop searching
						}

					// otherwise we're missing a leading entry in the main diagonal, which means det is zero
					// setting col to c makes it think there was a row of zeroes, and thus det of zero
					if (only_det) col = c;
				}
			found_entry:

				// if we didn't find a leading entry, the rest of the matrix is zeroes (so we're done)
				if (col == c)
				{
					if (det) *det = std::move(zero); // no leading entry means a row of all zero, which means a det of zero
					return row;
				}

				// make this row's leading entry a 1 via row division //

				// store lead entry for efficiency
				temp = operator()(row, col);
				// all elements to left of leading entry are zeros, so we can ignore them
				for (std::size_t i = col + 1; i < c; ++i) operator()(row, i) /= temp;
				// setting leading entry to 1 is more efficient and prevents rounding errors
				operator()(row, col) = (T)1;
				_det *= temp; // account for this entry in the determinant

				// use this row to eliminate the higher rows //

				// only do this if requested
				if (kill_upper)
				{
					for (std::size_t i = 0; i < row; ++i)
					{
						// store factor to subtract by for efficiency
						temp = operator()(i, col);
						if (temp != zero)
						{
							// all source elements to left of col are zero, so we can ignore them
							for (std::size_t j = col + 1; j < c; ++j) operator()(i, j) -= operator()(row, j) * temp;
							// setting entry to 0 is more efficient and prevents rounding errors
							operator()(i, col) = zero;
						}
					}
				}

				// use this row to eliminate the lower rows //

				for (std::size_t i = row + 1; i < r; ++i)
				{
					// store factor to subtract by for efficiency
					temp = operator()(i, col);
					if (temp != zero)
					{
						// all source elements to left of col are zero, so we can ignore them
						for (std::size_t j = col + 1; j < c; ++j) operator()(i, j) -= operator()(row, j) * temp;
						// setting entry to 0 is more efficient and prevents rounding errors
						operator()(i, col) = zero;
					}
				}
			}

			// export the calculated determinant
			if (det) *det = std::move(_det);
			return r;
		}

	public: // -- element indexing -- //

		// gets the element at the specified row and col
		T &operator()(std::size_t row, std::size_t col) { return data[row * c + col]; }
		const T &operator()(std::size_t row, std::size_t col) const { return data[row * c + col]; }

		// gets the element at the specified row and col, with bounds checking
		T &at(std::size_t row, std::size_t col)
		{
			if (row >= r || col >= c) throw std::out_of_range("matrix position out of bounds");
			return operator()(row, col);
		}
		const T &at(std::size_t row, std::size_t col) const
		{
			if (row >= r || col >= c) throw std::out_of_range("matrix position out of bounds");
			return operator()(row, col);
		}

	public: // -- row indexing -- //

		// gets a view into the specified row
		row_view operator[](std::size_t row) { return { data.data() + row * c, c }; }
		const_row_view operator[](std::size_t row) const { return { data.data() + row * c, c }; }

		// gets a view into the specified row, with bounds checking
		row_view at(std::size_t row)
		{
			if (row >= r) throw std::out_of_range("row out of bounds");
			return operator[](row);
		}
		const_row_view at(std::size_t row) const
		{
			if (row >= r) throw std::out_of_range("row out of bounds");
			return operator[](row);
		}

	public: // -- iterators -- //

		// iterates through each row in the matrix
		row_iter begin() noexcept { return { data.data(), c }; }
		row_iter end() noexcept { return { data.data() + data.size(), c }; }

		const_row_iter begin() const noexcept { return { data.data(), c }; }
		const_row_iter end() const noexcept { return { data.data() + data.size(), c }; }

		const_row_iter cbegin() const noexcept { return begin(); }
		const_row_iter cend() const noexcept { return end(); }

		// iteratres through each row in the matrix in reverse order
		reverse_row_iter rbegin() noexcept { return { end() }; }
		reverse_row_iter rend() noexcept { return { begin() }; }

		const_reverse_row_iter rbegin() const noexcept { return { end() }; }
		const_reverse_row_iter rend() const noexcept { return { begin() }; }

		const_reverse_row_iter crbegin() const noexcept { return rbegin(); }
		const_reverse_row_iter crend() const noexcept { return rend(); }

	public: // -- views -- //

		// views the entire matrix as a flattened array (row major order)
		flat_view flat() noexcept { return { data.data(), data.size() }; }
		const_flat_view flat() const noexcept { return { data.data(), data.size() }; }

	public: // -- ctor / dtor / asgn -- //

		// creates an empty matrix (0x0)
		matrix() noexcept = default;

		matrix(const matrix &other) = default;
		matrix &operator=(const matrix &other) = default;

		// constructs from another matrix by using its allocated resources.
		// other is guaranteed to be empty() after this operation.
		matrix(matrix &&other) noexcept : data(std::move(other.data)), r(std::exchange(other.r, 0)), c(std::exchange(other.c, 0)) {}
		// copies from another matrix via move semantics.
		// other is guaranteed to be empty() after this operation.
		matrix &operator=(matrix &&other) noexcept(noexcept(*(decltype(data)*)0 = std::declval<decltype(data)>()))
		{
			if (this != &other)
			{
				data = std::move(other.data);
				r = other.r;
				c = other.c;

				other.clear();
			}
			return *this;
		}

		// creates a matrix with the specified dimensions.
		// elements are copied from val.
		matrix(std::size_t rows, std::size_t cols, const T &val)
		{
			resize_flat(rows, cols, val);
		}
		// creates a matrix with the specified dimensions.
		// elements are defaulted to zero.
		matrix(std::size_t rows, std::size_t cols) : matrix(rows, cols, (T)0) {}

		// creates a matrix with the given values
		matrix(std::initializer_list<std::initializer_list<T>> vals)
		{
			*this = vals;
		}
		// assigns the given values to this matrix
		matrix &operator=(std::initializer_list<std::initializer_list<T>> vals)
		{
			clear();

			const std::size_t rows = vals.size();
			if (rows == 0) return *this;

			auto row = vals.begin();
			const std::size_t cols = row->size();

			data.reserve(rows * cols);
			for (; row != vals.end(); ++row)
			{
				if (row->size() != cols) throw std::invalid_argument("matrix initializer list was not rectangular");
				for (const T &val : *row) data.push_back(val);
			}
			if (!data.empty())
			{
				r = rows;
				c = cols;
			}

			return *this;
		}

	public: // -- shape utilities -- //

		// gets the number of rows in the matrix
		std::size_t rows() const { return r; }
		// gets the number of columns in the matrix
		std::size_t cols() const { return c; }

		// returns the total number of elements in the matrix (rows * cols)
		std::size_t size() const { return data.size(); }
		// returns the current capacity of the matrix (as a flattened array)
		std::size_t capacity() const { return data.capacity(); }

		// requests the matrix to set aside space for the specified number of total elements.
		// contents of the matrix are preserved.
		void reserve(std::size_t count)
		{
			data.reserve(count);
		}
		// requests the matrix to remove unused allocated space (non-binding).
		// contents of the matrix are preserved.
		void shrink_to_fit()
		{
			data.shrink_to_fit();
		}

		// returns true if the matrix is empty (0x0)
		bool empty() const { return data.empty(); }
		// sets the matrix to the empty state (0x0)
		void clear()
		{
			data.clear();
			r = 0;
			c = 0;
		}

		// resizes the matrix, preserving the flattened positions of any pre-existing elements present in both sizes.
		// shrinking truncates elements. expanding adds copies of val.
		// if either new dimension is zero this is equivalent to clear().
		void resize_flat(std::size_t rows, std::size_t cols, const T &val)
		{
			if (rows == 0 || cols == 0) { clear(); return; }
			data.resize(rows * cols, val);
			r = rows;
			c = cols;
		}
		// as resize_flat() but uses (T)0 as the fill value
		void resize_flat(std::size_t rows, std::size_t cols)
		{
			resize_flat(rows, cols, (T)0);
		}

		// resizes the number of columns in this matrix, preserving the (row, col) positions of any pre-existing elements present in both sizes.
		// shrinking truncates elements from the right. expanding adds copies of val to the right.
		// if cols is zero this is equivalent to clear().
		void resize_cols(std::size_t cols, const T &val)
		{
			if (cols == 0) clear();
			else if (cols > c)
			{
				data.resize(r * cols); // expand the array to needed size
				for (std::size_t i = r; i-- > 1; ) // for each row in reverse order (except the first)
				{
					T *const src = data.data() + i * c;
					T *const dest = data.data() + i * cols;

					for (std::size_t j = cols; j-- > c; ) dest[j] = val; // pad the right with copies of val
					for (std::size_t j = c; j-- > 0; ) dest[j] = std::move(src[j]); // shuffle elements up to their correct new flattened positions (backwards)
				}
				for (std::size_t j = cols; j-- > c; ) data[j] = val; // finish padding first row (separate from above loop to avoid glorified no-op element moves when dest == src)
				c = cols; // update number of columns
			}
			else if (cols < c)
			{
				for (std::size_t i = 1; i < r; ++i) // for each row (except the first)
				{
					T *const src = data.data() + i * c;
					T *const dest = data.data() + i * cols;

					std::move(src, src + cols, dest); // shuffle elements down to their correct new flattened positions
				}
				data.resize(r * cols); // shrink down the array
				c = cols; // update number of columns
			}
		}
		// as resize_cols() but uses (T)0 as the fill value
		void resize_cols(std::size_t cols)
		{
			resize_cols(cols, (T)0);
		}

		// resizes the number of rows in this matrx, preserving the (row, col) positions of an pre-existing elements present in both sizes.
		// shrinking truncates elements from the bottom. expanding adds copies of val to the bottom.
		// if rows is zero this is equivalent to clear().
		void resize_rows(std::size_t rows, const T &val)
		{
			if (rows == 0) clear();
			else
			{
				data.resize(rows * c, val);
				r = rows;
			}
		}
		// as resize_rows() but uses (T)0 as the fill value
		void resize_rows(std::size_t rows)
		{
			resize_rows(rows, (T)0);
		}

		// resizes the matrix to the given dimensions, preserving the (row, col) positions of any pre-existing elements present in both sizes.
		// the result is as if by calling both resize_rows() and resize_cols().
		// if either new dimension is zero this is equivalent to clear().
		void resize(std::size_t rows, std::size_t cols, const T &val)
		{
			if (rows == 0 || cols == 0) clear(); // if either is zero just clear matrix to avoid unnecessary work from naive row/col resizing
			else if (rows < r) // otherwise resize rows/cols in such a way that resize_cols has to deal with as few rows as possible
			{
				resize_rows(rows, val);
				resize_cols(cols, val);
			}
			else
			{
				resize_cols(cols, val);
				resize_rows(rows, val);
			}
		}
		// as resize() but uses (T)0 as the fill value
		void resize(std::size_t rows, std::size_t cols)
		{
			resize(rows, cols, (T)0);
		}

		// appends the columns of another matrix to the right of this one (allows self-appending).
		// throws matrix_size_error if the two matrices differ in number of rows.
		// if other is passed as rvalue, other is guaranteed to be empty() after this operation.
		// if self-appending, other should not be passed as rvalue.
		template<typename U, std::enable_if_t<std::is_same_v<std::remove_cv_t<std::remove_reference_t<U>>, matrix>, int> = 0>
		void append_cols(U &&other)
		{
			if (r != other.r) throw matrix_size_error("attempt to append_cols on matrices of different row counts");
			const std::size_t cols = c + other.c; // number of resulting columns

			constexpr bool is_rvalue = std::is_same_v<U, matrix> || std::is_same_v<U, matrix&&>;
			typedef std::conditional_t<is_rvalue, mover, noop> migrator;

			data.resize(r * cols); // expand the array to needed size
			for (std::size_t i = r; i-- > 1; ) // for each row in reverse order (except the first)
			{
				T *const src = data.data() + i * c;
				T *const dest = data.data() + i * cols;

				for (std::size_t j = c; j-- > 0; ) dest[j] = std::move(src[j]); // shuffle elements up to their correct new flattened positions (backwards)
			}
			const std::size_t other_c_fix = this != &other ? other.c : cols; // account for self-append column shifting in right fill
			for (std::size_t i = 0; i < r; ++i) // for each row 
			{
				auto *const src = other.data.data() + i * other_c_fix;
				T *const dest = data.data() + i * cols;

				for (std::size_t j = c; j < cols; ++j) dest[j] = migrator{}(src[j - c]); // fill right side with copies/moves of other's elements
			}
			c = cols; // update number of columns

			// if this is an rvalue append operation, we have special conditions
			if constexpr (is_rvalue)
			{
				assert(this != &other && "do not pass rvalue for self-appending");
				other.clear();
			}
		}





















	

		// computes the tensor product of the left matrix by the right matrix
		static matrix tensor_product(const matrix &left, const matrix &right)
		{
			matrix res(left.r * right.r, left.c * right.c);

			// compute the tensor product entries
			for (std::size_t left_row = 0; left_row < left.r; ++left_row)
				for (std::size_t left_col = 0; left_col < left.c; ++left_col)
				{
					// hold on to row position in result (shorthand)
					const std::size_t row = left_row * right.r;
					const std::size_t col = left_col * right.c;

					for (std::size_t right_row = 0; right_row < right.r; ++right_row)
						for (std::size_t right_col = 0; right_col < right.c; ++right_col)
							res(row + right_row, col + right_col) = left(left_row, left_col) * right(right_row, right_col);
				}

			return res;
		}

	public: // -- elementary row operations -- //

		// swaps rows a and b
		void swapRows(std::size_t a, std::size_t b)
		{
			using std::swap;
			if (a == b) return;

			for (std::size_t i = 0; i < c; ++i)
				swap(operator()(a, i), operator()(b, i));
		}

	public: // -- operations -- //

		// returns the determinant of the matrix.
		// throws matrix_size_error if matrix is non-square.
		T det() &&
		{
			if (r != c) throw matrix_size_error("attempt to take determinant of non-square matrix");
			if (r == 0) return (T)1; // determinant of 0x0 matrix is 1 to appease various other general properties

			T res;
			reduce(false, &res, true);
			return res;
		}
		T det() const&
		{
			if (r != c) throw matrix_size_error("attempt to take determinant of non-square matrix");
			if (r == 0) return (T)1; // determinant of 0x0 matrix is 1 to appease various other general properties

			matrix cpy = *this;
			T res;
			cpy.reduce(false, &res, true);
			return res;
		}

		// returns true if the matrix is invertible
		bool invertible() const& { return r == c && det() != (T)0; }
		bool invertible() && { return r == c && std::move(*this).det() != (T)0; }

		// returns the <row><col> submatrix (i.e. the result of removing row <row> and col <col>)
		// throws matrix_size_error if matrix is empty
		matrix submatrix(std::size_t row, std::size_t col) const
		{
			if (r == 0) throw matrix_size_error("Cannot take a submatrix from an empty matrix");

			// allocate the result
			matrix res(r - 1, c - 1);

			// fill its entries
			for (std::size_t _row = 0; _row < row; ++_row)
			{
				for (std::size_t _col = 0; _col < col; ++_col)
					res(_row, _col) = operator()(_row, _col);
				for (std::size_t _col = col + 1; _col < c; ++_col)
					res(_row, _col - 1) = operator()(_row, _col);
			}
			for (std::size_t _row = row + 1; _row < r; ++_row)
			{
				for (std::size_t _col = 0; _col < col; ++_col)
					res(_row - 1, _col) = operator()(_row, _col);
				for (std::size_t _col = col + 1; _col < c; ++_col)
					res(_row - 1, _col - 1) = operator()(_row, _col);
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
			return ((row + col) & 1 ? -1 : 1) * at(row, col) * minor(row, col);
		}

		// returns the cofactor matrix of this matrix
		// throws matrix_size_error if non-square matrix or empty
		// also throws by cofactor()
		// complexity: O(n^2) if invertible, otherwise O(n^4)
		matrix cofactormatrix() const
		{
			if (r != c) throw matrix_size_error("Cannot take cofactor of non-square matrix");
			if (r == 0) throw matrix_size_error("Cannot take cofactor of empty matrix");

			// allocate space for result
			matrix res(r, r);
			T _det;

			// definition of adjugate: A * adj(A) = det(A) * I
			// if A is invertible:     adj(A) = det(A) * inv(A)
			// then, since the adjugate is transpose of the cofactor matrix, just transpose adjugate
			if (inverse(res, &_det))
			{
				res *= _det;
				res.transpose();
			}
			// otherwise we must resort to the other definition: adj(A) = cofactormatrix(A).transpose()
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
		// throws matrix_size_error if non-square matrix or empty
		// also throws by cofactor()
		// complexity: O(n^2) if invertible, otherwise O(n^4)
		matrix adjugate() const
		{
			if (r != c) throw matrix_size_error("Cannot take adjugate of non-square matrix");
			if (r == 0) throw matrix_size_error("Cannot take adjugate of empty matrix");

			// allocate space for result
			matrix res(r, r);
			T _det;

			// definition of adjugate: A * adj(A) = det(A) * I
			// if A is invertible:     adj(A) = det(A) * inv(A)
			if (inverse(res, &_det)) res *= _det;
			// otherwise we must resort to the other definition: adj(A) = cofactormatrix(A).transpose()
			else
			{
				// fill with cofactors
				for (std::size_t row = 0; row < r; ++row)
					for (std::size_t col = 0; col < c; ++col)
						res(col, row) = cofactor(row, col);
			}

			return res;
		}

		// transposes the matrix
		matrix &transpose()
		{
			using std::swap;

			// as a special case, row/column vectors will have identical flattened structures
			if (r == 1 || c == 1) swap(r, c);
			// square matrices can use simple swaps
			else if (r == c)
			{
				for (std::size_t row = 0; row < r; ++row)
					for (std::size_t col = row + 1; col < c; ++col)
						swap(operator()(row, col), operator()(col, row));
			}
			// otherwise it's complicated. just copy a transposed version
			else *this = getTranspose();

			return *this;
		}
		// returns a copy of this matrix that has been transposed
		matrix getTranspose() const
		{
			// allocate the matrix
			matrix res(c, r);

			// copy the transposed entries
			for (std::size_t row = 0; row < r; ++row)
				for (std::size_t col = 0; col < c; ++col)
					res(col, row) = operator()(row, col);

			return res;
		}

		// puts the matrix into row echelon form. returns the rank of the matrix
		// optionally also returns the determinant of the largest square matrix that this matrix can contain
		// complexity: O(nm)
		std::size_t REF(T *det = nullptr) { return reduce(false, det, false); }
		// puts the matrix into reduced row echelon form. returns the rank of the matrix
		// complexity: O(nm)
		std::size_t RREF(T *det = nullptr) { return reduce(true, det, false); }

		// attempts to find the inverse of the matrix. if successful, stores inverse in <dest> and returns true
		// optionally also returns the determinant of the largest square matrix that this matrix can contain
		// throws matrix_size_error if matrix is not square or is empty
		// complexity: O(n^2)
		bool inverse(matrix &dest, T *det = nullptr) const
		{
			if (r != c) throw matrix_size_error("Only square matrices are invertible");
			if (r == 0) throw matrix_size_error("Cannot take the inverse of an empty matrix");

			matrix util(r, 2 * c); // make a wide matrix

			// fill left side with current matrix and right side with identity matrix
			for (std::size_t row = 0; row < r; ++row)
				for (std::size_t col = 0; col < c; ++col)
				{
					util(row, col) = (*this)(row, col);
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
		bool invert(T *det = nullptr) { return inverse(*this, det); }
	};

	// -- operator definitions -- //

	template<typename T, typename Allocator>
	matrix<T, Allocator> &operator+=(matrix<T, Allocator> &a, const matrix<T, Allocator> &b)
	{
		if (a.rows() != b.rows() || a.cols() != b.cols()) throw matrix_size_error("matrix addition requires the matrices be of same size");

		std::for_each_n(exec_policy<T>(), val_iter<std::size_t>{0}, a.size(),
			[f1 = a.flat(), f2 = b.flat()](std::size_t i) { return f1[i] += f2[i]; });

		return a;
	}

	template<typename T, typename Allocator> matrix<T, Allocator> operator+(const matrix<T, Allocator> &a, const matrix<T, Allocator> &b) { auto res = a; res += b; return res; }
	template<typename T, typename Allocator> matrix<T, Allocator> operator+(matrix<T, Allocator> &&a, const matrix<T, Allocator> &b) { a += b; return std::move(a); }
	template<typename T, typename Allocator> matrix<T, Allocator> operator+(const matrix<T, Allocator> &a, matrix<T, Allocator> &&b) { b += a; return std::move(b); }
	template<typename T, typename Allocator> matrix<T, Allocator> operator+(matrix<T, Allocator> &&a, matrix<T, Allocator> &&b) { a += b; return std::move(a); }

	template<typename T, typename Allocator>
	matrix<T, Allocator> &operator-=(matrix<T, Allocator> &a, const matrix<T, Allocator> &b)
	{
		if (a.rows() != b.rows() || a.cols() != b.cols()) throw matrix_size_error("matrix subtraction requires the matrices be of same size");

		std::for_each_n(exec_policy<T>(), val_iter<std::size_t>{0}, a.size(), [f1 = a.flat(), f2 = b.flat()](std::size_t i) { return f1[i] -= f2[i]; });

		return a;
	}

	template<typename T, typename Allocator> matrix<T, Allocator> operator-(const matrix<T, Allocator> &a, const matrix<T, Allocator> &b) { auto res = a; res -= b; return res; }
	template<typename T, typename Allocator> matrix<T, Allocator> operator-(matrix<T, Allocator> &&a, const matrix<T, Allocator> &b) { a -= b; return std::move(a); }

	template<typename T, typename Allocator> matrix<T, Allocator> &operator*=(matrix<T, Allocator> &m, const T &r)
	{
		for (T &v : m.flat()) v *= r;
		return m;
	}

	template<typename T, typename Allocator> matrix<T, Allocator> operator*(const matrix<T, Allocator> &m, const T &f) { auto res = m; res *= f; return res; }
	template<typename T, typename Allocator> matrix<T, Allocator> operator*(const T &f, const matrix<T, Allocator> &m) { auto res = m; res *= f; return res; }
	template<typename T, typename Allocator> matrix<T, Allocator> operator*(matrix<T, Allocator> &&m, const T &f) { m *= f; return std::move(m); }
	template<typename T, typename Allocator> matrix<T, Allocator> operator*(const T &f, matrix<T, Allocator> &&m) { m *= f; return std::move(m); }

	template<typename T, typename Allocator>
	matrix<T, Allocator> &operator/=(matrix<T, Allocator> &m, T f)
	{
		for (std::size_t row = 0; row < m.rows(); ++row)
			for (std::size_t col = 0; col < m.cols(); ++col)
				m(row, col) /= f;

		return m;
	}

	template<typename T, typename Allocator> matrix<T, Allocator> operator/(const matrix<T, Allocator> &m, T f) { auto res = m; res /= f; return res; }
	template<typename T, typename Allocator> matrix<T, Allocator> operator/(matrix<T, Allocator> &&m, T f) { m /= f; return std::move(m); }

	template<typename T, typename Allocator>
	matrix<T, Allocator> operator*(const matrix<T, Allocator> &lhs, const matrix<T, Allocator> &rhs)
	{
		if (lhs.cols() != rhs.rows()) throw matrix_size_error("matrix multiplication requires lhs #cols equal rhs #rows");

		matrix<T, Allocator> res(lhs.rows(), rhs.cols());

		std::for_each_n(exec_policy<T>(), val_iter<std::size_t>{0}, res.size(), [&res, &lhs, &rhs](std::size_t idx)
		{
			const std::size_t row = idx / res.cols();
			const std::size_t col = idx % res.cols();

			T dot{ 0 };
			for (std::size_t i = 0; i < lhs.cols(); ++i) dot += lhs(row, i) * rhs(i, col);
			res.flat()[idx] = std::move(dot);
		});

		return res;
	}
	template<typename T, typename Allocator> matrix<T, Allocator> &operator*=(matrix<T, Allocator> &lhs, const matrix<T, Allocator> &rhs) { return lhs = lhs * rhs; }

	// -- cmp operator definitions -- //

	template<typename T, typename Allocator>
	bool operator==(const matrix<T, Allocator> &a, const matrix<T, Allocator> &b)
	{
		if (a.rows() != b.rows() || a.cols() != b.cols()) return false;

		return std::all_of(exec_policy<T>(), val_iter<std::size_t>{0}, val_iter<std::size_t>{a.size()}, [f1 = a.flat(), f2 = b.flat()](std::size_t i) { return f1[i] == f2[i]; });
	}
	template<typename T, typename Allocator>
	bool operator!=(const matrix<T, Allocator> &a, const matrix<T, Allocator> &b) { return !(a == b); }
}

template<typename T, typename Allocator = std::allocator<T>>
using matrix = matrix_impl::matrix<T, Allocator>;

using matrix_impl::matrix_size_error;

#endif
