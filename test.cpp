#include <iostream>
#include <iomanip>
#include <chrono>
#include <random>
#include <functional>
#include <type_traits>
#include <iterator>

#undef NDEBUG
#include <cassert>

#include "matrix.h"

// -------------------------------------------------------------------------------------------------------------------------------------------

#define assert_throws(expr, ex) { try { (void)(expr); assert(false && "did not throw"); } catch (const ex &) {} catch (...) { assert(false && "threw wrong type"); } }

// -------------------------------------------------------------------------------------------------------------------------------------------

static_assert(std::is_same_v<std::iterator_traits<matrix_impl::val_iter<std::size_t>>::iterator_category, std::random_access_iterator_tag>, "val iterator is not random access");
static_assert(std::is_same_v<std::iterator_traits<matrix<double>::iterator>::iterator_category, std::random_access_iterator_tag>, "matrix iterator is not random access");

static_assert(std::is_same_v<std::iterator_traits<matrix<double>::row_view::iterator>::iterator_category, std::random_access_iterator_tag>, "row view iterator is not random access");
static_assert(std::is_same_v<std::iterator_traits<matrix<double>::const_row_view::iterator>::iterator_category, std::random_access_iterator_tag>, "row view iterator is not random access");

static_assert(std::is_same_v<std::iterator_traits<matrix<double>::iterator>::reference, matrix<double>::row_view>, "matrix iterator is not a row iterator");

// -------------------------------------------------------------------------------------------------------------------------------------------

static_assert(std::is_trivially_default_constructible_v<matrix<double>::flat_view>, "matrix flat views should be trivial");
static_assert(std::is_trivially_default_constructible_v<matrix<double>::const_flat_view>, "matrix flat views should be trivial");

static_assert(std::is_trivially_destructible_v<matrix<double>::flat_view>, "matrix flat views should be trivial");
static_assert(std::is_trivially_destructible_v<matrix<double>::const_flat_view>, "matrix flat views should be trivial");

static_assert(std::is_trivially_copy_constructible_v<matrix<double>::flat_view>, "matrix flat views should be trivial");
static_assert(std::is_trivially_copy_constructible_v<matrix<double>::const_flat_view>, "matrix flat views should be trivial");

static_assert(std::is_trivially_copy_assignable_v<matrix<double>::flat_view>, "matrix flat views should be trivial");
static_assert(std::is_trivially_copy_assignable_v<matrix<double>::const_flat_view>, "matrix flat views should be trivial");

static_assert(std::is_convertible_v<matrix<double>::flat_view, matrix<double>::const_flat_view>, "matrix flat views should allow standard cv conversions");
static_assert(!std::is_convertible_v<matrix<double>::const_flat_view, matrix<double>::flat_view>, "matrix flat views should allow standard cv conversions");

static_assert(!std::is_convertible_v<matrix<double>::const_flat_view, matrix<double>::flat_view>, "matrix flat views should allow standard cv conversions");
static_assert(std::is_convertible_v<matrix<double>::flat_view, matrix<double>::const_flat_view>, "matrix flat views should allow standard cv conversions");

// -------------------------------------------------------------------------------------------------------------------------------------------

static_assert(std::is_trivially_default_constructible_v<matrix<double>::row_view>, "matrix row views should be trivial");
static_assert(std::is_trivially_default_constructible_v<matrix<double>::const_row_view>, "matrix row views should be trivial");

static_assert(std::is_trivially_destructible_v<matrix<double>::row_view>, "matrix row views should be trivial");
static_assert(std::is_trivially_destructible_v<matrix<double>::const_row_view>, "matrix row views should be trivial");

static_assert(std::is_trivially_copy_constructible_v<matrix<double>::row_view>, "matrix row views should be trivial");
static_assert(std::is_trivially_copy_constructible_v<matrix<double>::const_row_view>, "matrix row views should be trivial");

static_assert(std::is_trivially_copy_assignable_v<matrix<double>::row_view>, "matrix row views should be trivial");
static_assert(std::is_trivially_copy_assignable_v<matrix<double>::const_row_view>, "matrix row views should be trivial");

static_assert(std::is_convertible_v<matrix<double>::row_view, matrix<double>::const_row_view>, "matrix row views should allow standard cv conversions");
static_assert(!std::is_convertible_v<matrix<double>::const_row_view, matrix<double>::row_view>, "matrix row views should allow standard cv conversions");

static_assert(!std::is_convertible_v<matrix<double>::const_row_view, matrix<double>::row_view>, "matrix row views should allow standard cv conversions");
static_assert(std::is_convertible_v<matrix<double>::row_view, matrix<double>::const_row_view>, "matrix row views should allow standard cv conversions");

// -------------------------------------------------------------------------------------------------------------------------------------------

static_assert(std::is_trivially_default_constructible_v<matrix<double>::row_iter>, "matrix row iter should be trivial");
static_assert(std::is_trivially_default_constructible_v<matrix<double>::const_row_iter>, "matrix row iter should be trivial");

static_assert(std::is_trivially_destructible_v<matrix<double>::row_iter>, "matrix row iter should be trivial");
static_assert(std::is_trivially_destructible_v<matrix<double>::const_row_iter>, "matrix row iter should be trivial");

static_assert(std::is_trivially_copy_constructible_v<matrix<double>::row_iter>, "matrix row iter should be trivial");
static_assert(std::is_trivially_copy_constructible_v<matrix<double>::const_row_iter>, "matrix row iter should be trivial");

static_assert(std::is_trivially_copy_assignable_v<matrix<double>::row_iter>, "matrix row iter should be trivial");
static_assert(std::is_trivially_copy_assignable_v<matrix<double>::const_row_iter>, "matrix row iter should be trivial");

static_assert(std::is_convertible_v<matrix<double>::row_iter, matrix<double>::const_row_iter>, "matrix row iter should allow standard cv conversions");
static_assert(!std::is_convertible_v<matrix<double>::const_row_iter, matrix<double>::row_iter>, "matrix row iter should allow standard cv conversions");

static_assert(!std::is_convertible_v<matrix<double>::const_row_iter, matrix<double>::row_iter>, "matrix row iter should allow standard cv conversions");
static_assert(std::is_convertible_v<matrix<double>::row_iter, matrix<double>::const_row_iter>, "matrix row iter should allow standard cv conversions");

// -------------------------------------------------------------------------------------------------------------------------------------------

static_assert(std::is_trivially_default_constructible_v<matrix<double>::row_view::iterator>, "row view iterator should be trivial");
static_assert(std::is_trivially_default_constructible_v<matrix<double>::const_row_view::iterator>, "row view iterator should be trivial");

static_assert(std::is_trivially_destructible_v<matrix<double>::row_view::iterator>, "row view iterator should be trivial");
static_assert(std::is_trivially_destructible_v<matrix<double>::const_row_view::iterator>, "row view iterator should be trivial");

static_assert(std::is_trivially_copy_constructible_v<matrix<double>::row_view::iterator>, "row view iterator should be trivial");
static_assert(std::is_trivially_copy_constructible_v<matrix<double>::const_row_view::iterator>, "row view iterator should be trivial");

static_assert(std::is_trivially_copy_assignable_v<matrix<double>::row_view::iterator>, "row view iterator should be trivial");
static_assert(std::is_trivially_copy_assignable_v<matrix<double>::const_row_view::iterator>, "row view iterator should be trivial");

// -------------------------------------------------------------------------------------------------------------------------------------------

static_assert(std::is_same_v<decltype(*(*(matrix<double>*)0).begin()), matrix<double>::row_view>, "row iter does not return row view");
static_assert(std::is_same_v<decltype(*(*(const matrix<double>*)0).begin()), matrix<double>::const_row_view>, "row iter does not return row view");
static_assert(std::is_same_v<decltype(*(*(matrix<double>*)0).cbegin()), matrix<double>::const_row_view>, "row iter does not return row view");
static_assert(std::is_same_v<decltype(*(*(const matrix<double>*)0).cbegin()), matrix<double>::const_row_view>, "row iter does not return row view");

static_assert(std::is_same_v<decltype(*(*(matrix<double>*)0).rbegin()), matrix<double>::row_view>, "rev row iter does not return row view");
static_assert(std::is_same_v<decltype(*(*(const matrix<double>*)0).rbegin()), matrix<double>::const_row_view>, "rev row iter does not return row view");
static_assert(std::is_same_v<decltype(*(*(matrix<double>*)0).crbegin()), matrix<double>::const_row_view>, "rev row iter does not return row view");
static_assert(std::is_same_v<decltype(*(*(const matrix<double>*)0).crbegin()), matrix<double>::const_row_view>, "rrevow iter does not return row view");

// -------------------------------------------------------------------------------------------------------------------------------------------

static_assert(std::is_same_v<decltype(std::declval<matrix<double>::row_view>()[0]), double&>, "row view does not give references");
static_assert(std::is_same_v<decltype(std::declval<matrix<double>::const_row_view>()[0]), const double&>, "const row view does not give references");

static_assert(std::is_same_v<decltype(std::declval<matrix<double>::row_view>().at(0)), double&>, "row view does not give references");
static_assert(std::is_same_v<decltype(std::declval<matrix<double>::const_row_view>().at(0)), const double&>, "const row view does not give references");

// -------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
struct intrinsic
{
	static constexpr T move_ctor = (T)5346;
	static constexpr T move_assign = (T)72853;

	static_assert(std::is_integral_v<T>, "T must be an integral type");

	T val;
	intrinsic() = default;
	explicit intrinsic(T v) : val(v) {}

	intrinsic(const intrinsic&) = default;
	intrinsic &operator=(const intrinsic&) = default;

	intrinsic(intrinsic &&other) noexcept : val(std::exchange(other.val, move_ctor)) {}
	intrinsic &operator=(intrinsic &&other) noexcept { val = std::exchange(other.val, move_assign); return *this; }

	operator T() const { return val; }

	friend intrinsic operator+(intrinsic a) { return a; }
	friend intrinsic operator-(intrinsic a) { return intrinsic{ a.val }; }

	friend intrinsic &operator+=(intrinsic &a, intrinsic b) { a.val += b.val; return a; }
	friend intrinsic operator+(intrinsic a, intrinsic b) { return intrinsic{ a.val + b.val }; }

	friend intrinsic &operator-=(intrinsic &a, intrinsic b) { a.val -= b.val; return a; }
	friend intrinsic operator-(intrinsic a, intrinsic b) { return intrinsic{ a.val - b.val }; }

	friend intrinsic &operator*=(intrinsic &a, intrinsic b) { a.val *= b.val; return a; }
	friend intrinsic operator*(intrinsic a, intrinsic b) { return intrinsic{ a.val * b.val }; }

	friend intrinsic &operator/=(intrinsic &a, intrinsic b) { a.val /= b.val; return a; }
	friend intrinsic operator/(intrinsic a, intrinsic b) { return intrinsic{ a.val / b.val }; }

	friend bool operator==(intrinsic a, intrinsic b) { return a.val == b.val; }
	friend bool operator!=(intrinsic a, intrinsic b) { return a.val != b.val; }
	friend bool operator<(intrinsic a, intrinsic b) { return a.val < b.val; }
	friend bool operator<=(intrinsic a, intrinsic b) { return a.val <= b.val; }
	friend bool operator>(intrinsic a, intrinsic b) { return a.val > b.val; }
	friend bool operator>=(intrinsic a, intrinsic b) { return a.val >= b.val; }
};

// -------------------------------------------------------------------------------------------------------------------------------------------

// outputs the matrix to the stream
template<typename T>
std::ostream &operator<<(std::ostream &ostr, const matrix<T> &m)
{
	for (auto row : m)
	{
		for (const T &val : row) ostr << std::setw(10) << val << ' ';
		ostr << '\n';
	}
	return ostr;
}

template<typename T> auto to_signed(T t) noexcept { return static_cast<std::make_signed_t<T>>(t); }
template<typename T> auto to_unsigned(T t) noexcept { return static_cast<std::make_unsigned_t<T>>(t); }

// -------------------------------------------------------------------------------------------------------------------------------------------

static auto &rng()
{
	static std::mt19937 r{ std::random_device{}() };
	return r;
}

void rand_fill(matrix<double> &m)
{
	auto &r = rng();
	std::uniform_real_distribution<double> dist(-1, 1);
	for (double &v : m.flat()) v = dist(r);
}

// -------------------------------------------------------------------------------------------------------------------------------------------

volatile std::size_t optimizer_killer = 0;

template<typename T>
void do_something_with_it(const matrix<T> &m) { optimizer_killer += (std::size_t)m.flat()[0]; }

template<typename T, typename F>
void benchmark_binary(const char *name, std::size_t reps, std::pair<std::size_t, std::size_t> d1, std::pair<std::size_t, std::size_t> d2, F f)
{
	matrix<T> m1(d1.first, d1.second);
	matrix<T> m2(d2.first, d2.second);

	rand_fill(m1);
	rand_fill(m2);

	const auto start = std::chrono::high_resolution_clock::now();

	for (std::size_t i = 0; i < reps; ++i)
	{
		decltype(auto) r = f(m1, m2);
		do_something_with_it(r);
	}

	const auto stop = std::chrono::high_resolution_clock::now();
	const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

	std::cout << std::left << std::setw(20) << name << ' ' << std::right << std::setw(10) << std::setprecision(2) << std::fixed << ((double)elapsed / reps) << " us avg (" << elapsed << ')' << std::endl;
}

// -------------------------------------------------------------------------------------------------------------------------------------------

void meta_tests()
{
	{
		typedef intrinsic<int> t;

		t a = (t)6;
		t b = a;
		assert(a == 6);
		assert(b == 6);
		assert(a == b);

		b = (t)7;
		assert(b == 7);
		assert(a != b);

		a = b;
		assert(a == 7);
		assert(a == b);

		a = (t)8;
		b = std::move(a);
		assert(b == 8);
		assert(a == t::move_assign);

		t c = std::move(b);
		assert(c == 8);
		assert(b == t::move_ctor);
	}
}

void basic_tests()
{
	matrix<std::size_t> a(10, 11);
	assert(a.rows() == 10);
	assert(a.cols() == 11);
	assert(a.size() == 110);

	const matrix<std::size_t>::flat_view f = a.flat();
	assert(f.size() == 110);

	for (std::size_t i = 0; i < f.size(); ++i) f[i] = i;

	for (std::size_t r = 0; r < a.rows(); ++r)
		for (std::size_t c = 0; c < a.cols(); ++c)
		{
			assert(a(r, c) == r * a.cols() + c);
		}

	{
		std::size_t v = 12;
		for (matrix<std::size_t>::row_view row : a)
			for (std::size_t &val : row)
			{
				val = v++;
			}
	}

	for (std::size_t r = 0; r < a.rows(); ++r)
		for (std::size_t c = 0; c < a.cols(); ++c)
		{
			assert(a.at(r, c) == r * a.cols() + c + 12);
		}

	assert_throws(a.at(a.rows(), 0), std::out_of_range);
	assert_throws(a.at(0, a.cols()), std::out_of_range);
	assert_throws(a.at(a.rows(), a.cols()), std::out_of_range);

	matrix<std::size_t> cpy = a;
	assert(cpy.rows() == a.rows());
	assert(cpy.cols() == a.cols());
	assert(cpy == a);
	assert(!cpy.empty());

	cpy.clear();
	assert(cpy.rows() == 0);
	assert(cpy.cols() == 0);
	assert(cpy.size() == 0);
	assert(cpy.empty());

	cpy = a;
	assert(cpy.rows() == a.rows());
	assert(cpy.cols() == a.cols());
	assert(cpy == a);
	assert(!cpy.empty());

	matrix<std::size_t> cpy2 = std::move(cpy);
	assert(cpy2.rows() == a.rows());
	assert(cpy2.cols() == a.cols());
	assert(cpy2 == a);
	assert(!cpy2.empty());
	assert(cpy.rows() == 0);
	assert(cpy.cols() == 0);
	assert(cpy.size() == 0);
	assert(cpy.empty());
	assert(cpy != a);

	cpy = a;
	assert(cpy.rows() == a.rows());
	assert(cpy.cols() == a.cols());
	assert(cpy == a);
	assert(!cpy.empty());

	cpy2 = std::move(cpy);
	assert(cpy2.rows() == a.rows());
	assert(cpy2.cols() == a.cols());
	assert(cpy2 == a);
	assert(!cpy2.empty());
	assert(cpy.rows() == 0);
	assert(cpy.cols() == 0);
	assert(cpy.size() == 0);
	assert(cpy.empty());
	assert(cpy != a);
}
void init_list_tests()
{
	{
		matrix<int> m;
		assert(m.rows() == 0);
		assert(m.cols() == 0);
		assert(m.size() == 0);
		assert(m.empty());
	}
	{
		matrix<int> m{};
		assert(m.rows() == 0);
		assert(m.cols() == 0);
		assert(m.size() == 0);
		assert(m.empty());
	}
	{
		matrix<int> m = {};
		assert(m.rows() == 0);
		assert(m.cols() == 0);
		assert(m.size() == 0);
		assert(m.empty());
	}
	{
		matrix<int> m{ {} };
		assert(m.rows() == 0);
		assert(m.cols() == 0);
		assert(m.size() == 0);
		assert(m.empty());
	}
	{
		matrix<int> m = { {} };
		assert(m.rows() == 0);
		assert(m.cols() == 0);
		assert(m.size() == 0);
		assert(m.empty());
	}
	{
		matrix<int> m{ 2, 4 }; // this is actually the size constructor (not 2d init list) - this is why 1d init list is not supported
		assert(m.rows() == 2);
		assert(m.cols() == 4);
		assert(m.size() == 8);
		assert(!m.empty());

		for (auto v : m.flat())
		{
			assert(v == 0);
		}
	}
	{
		matrix<int> m = { 2, 4 }; // this is actually the size constructor (not 2d init list) - this is why 1d init list is not supported
		assert(m.rows() == 2);
		assert(m.cols() == 4);
		assert(m.size() == 8);
		assert(!m.empty());

		for (auto v : m.flat())
		{
			assert(v == 0);
		}
	}
	{
		matrix<int> m{ 2, 4, 420 }; // this is actually the size constructor (not 2d init list) - this is why 1d init list is not supported
		assert(m.rows() == 2);
		assert(m.cols() == 4);
		assert(m.size() == 8);
		assert(!m.empty());

		for (auto row : m)
			for (int v : row)
			{
				assert(v == 420);
			}
	}
	{
		matrix<int> m = { 2, 4, 420 }; // this is actually the size constructor (not 2d init list) - this is why 1d init list is not supported
		assert(m.rows() == 2);
		assert(m.cols() == 4);
		assert(m.size() == 8);
		assert(!m.empty());

		for (auto row : m)
			for (int v : row)
			{
				assert(v == 420);
			}
	}
	{
		matrix<int> m{ { 1, 2, 3 } };
		assert(m.rows() == 1);
		assert(m.cols() == 3);
		assert(m.size() == 3);
		assert(!m.empty());
	}
	{
		matrix<int> m = { { 1, 2, 3 } };
		assert(m.rows() == 1);
		assert(m.cols() == 3);
		assert(m.size() == 3);
		assert(!m.empty());
	}
	{
		matrix<int> m{ { 1, 2, 3 }, { 4, 5, 6 } };
		assert(m.rows() == 2);
		assert(m.cols() == 3);
		assert(m.size() == 6);
		assert(!m.empty());

		assert(m.at(0, 0) == 1);
		assert(m.at(0, 1) == 2);
		assert(m.at(0, 2) == 3);
		assert(m.at(1, 0) == 4);
		assert(m.at(1, 1) == 5);
		assert(m.at(1, 2) == 6);
		assert_throws(m.at(2, 0), std::out_of_range);
		assert_throws(m.at(0, 3), std::out_of_range);
		assert_throws(m.at(1, 3), std::out_of_range);
		assert_throws(m.at(2, 3), std::out_of_range);

		assert(m(0, 0) == 1);
		assert(m(0, 1) == 2);
		assert(m(0, 2) == 3);
		assert(m(1, 0) == 4);
		assert(m(1, 1) == 5);
		assert(m(1, 2) == 6);

		assert(m.at(0).size() == 3);
		assert(m.at(1).size() == 3);

		assert(m.at(0).cols() == 3);
		assert(m.at(1).cols() == 3);

		assert(m.at(0).at(0) == 1);
		assert(m.at(0).at(1) == 2);
		assert(m.at(0).at(2) == 3);
		assert(m.at(1).at(0) == 4);
		assert(m.at(1).at(1) == 5);
		assert(m.at(1).at(2) == 6);

		assert(m[0].at(0) == 1);
		assert(m[0].at(1) == 2);
		assert(m[0].at(2) == 3);
		assert(m[1].at(0) == 4);
		assert(m[1].at(1) == 5);
		assert(m[1].at(2) == 6);
		assert_throws(m[1].at(3), std::out_of_range);

		assert(m.at(0)[0] == 1);
		assert(m.at(0)[1] == 2);
		assert(m.at(0)[2] == 3);
		assert(m.at(1)[0] == 4);
		assert(m.at(1)[1] == 5);
		assert(m.at(1)[2] == 6);
		assert_throws(m.at(2)[0], std::out_of_range);

		assert(m[0][0] == 1);
		assert(m[0][1] == 2);
		assert(m[0][2] == 3);
		assert(m[1][0] == 4);
		assert(m[1][1] == 5);
		assert(m[1][2] == 6);
	}
	{
		matrix<int> m = { { 1, 2, 3 }, { 4, 5, 6 } };
		assert(m.rows() == 2);
		assert(m.cols() == 3);
		assert(m.size() == 6);
		assert(!m.empty());

		assert(m.at(0, 0) == 1);
		assert(m.at(0, 1) == 2);
		assert(m.at(0, 2) == 3);
		assert(m.at(1, 0) == 4);
		assert(m.at(1, 1) == 5);
		assert(m.at(1, 2) == 6);
	}
	{
		matrix<int> m{ { 1, 0, 0, 0 }, { 0, 1, 0, 0 }, { 0, 0, 1, 0 }, { 0, 0, 0, 1 } };
		assert(m.rows() == 4);
		assert(m.cols() == 4);
		assert(m.size() == 16);
		assert(!m.empty());

		for (std::size_t i = 0; i < m.rows(); ++i)
			for (std::size_t j = 0; j < m.cols(); ++j)
			{
				assert(m.at(i, j) == (i == j));
			}
	}
	{
		matrix<int> m = { { 1, 0, 0, 0 }, { 0, 1, 0, 0 }, { 0, 0, 1, 0 }, { 0, 0, 0, 1 } };
		assert(m.rows() == 4);
		assert(m.cols() == 4);
		assert(m.size() == 16);
		assert(!m.empty());

		for (std::size_t i = 0; i < m.rows(); ++i)
			for (std::size_t j = 0; j < m.cols(); ++j)
			{
				assert(m.at(i, j) == (i == j));
			}
	}
	{
		assert_throws((matrix<int>{ {1, 2, 3}, { 4, 5 } }), std::invalid_argument);
		assert_throws((matrix<int>{ {1, 2, 3}, { 4, 5, 6, 7 } }), std::invalid_argument);
		assert_throws((matrix<int>{ {1, 2, 3}, {} }), std::invalid_argument);
		assert_throws((matrix<int>{ { }, { 1, 2, 3 } }), std::invalid_argument);
	}
}
void iter_tests()
{
	static_assert(std::is_same_v<std::iterator_traits<matrix<int>::iterator>::iterator_category, std::random_access_iterator_tag>, "matrix iter is not random access");
	static_assert(std::is_same_v<std::iterator_traits<matrix<int>::const_iterator>::iterator_category, std::random_access_iterator_tag>, "matrix const iter is not random access");

	matrix<int> m(18, 22);
	assert(m.rows() == 18);
	assert(m.cols() == 22);
	assert(!m.empty());

	assert(m.begin() == m.begin());
	assert((m.begin() + 0) == m.begin());
	assert((m.begin() - 0) == m.begin());
	assert((m.begin() + 1) != m.begin());
	for (std::size_t i = 1; i < m.rows(); ++i)
	{
		auto v = m.begin() + i;

		assert(v != m.begin());
		assert(v == v);
		assert(v != m.end());

		assert((v + 0) == v);
		assert((v - 0) == v);
		assert((v + 1) != v);
		assert((v - 1) != v);

		assert((std::size_t)(v - m.begin()) == i);
		assert((std::size_t)((v + 0) - m.begin()) == i);
		assert((std::size_t)((v - 0) - m.begin()) == i);
		assert((std::size_t)((v + 1) - m.begin()) == i + 1);
		assert((std::size_t)((v - 1) - m.begin()) == i - 1);

		assert((std::size_t)std::distance(m.begin(), v) == i);
		assert((std::size_t)std::distance(m.begin(), v + 0) == i);
		assert((std::size_t)std::distance(m.begin(), v - 0) == i);
		assert((std::size_t)std::distance(m.begin(), v + 1) == i + 1);
		assert((std::size_t)std::distance(m.begin(), v - 1) == i - 1);

		auto v2 = v;
		assert(v2 == v);
		auto v3 = v2++;
		assert(v3 == v);
		assert(v2 != v);
		assert(v2 != v3);

		auto v4 = v2--;
		assert(v2 == v);
		assert(v2 == v3);
		assert(v4 != v2);
		assert((v4 - v) == 1);
		assert(std::distance(v, v4) == 1);

		auto v5 = --v2;
		assert(v5 == v2);
		assert(v5 != v);
		assert(v2 != v);
		assert((v5 - v) == -1);
		assert(std::distance(v, v5) == -1);
		assert(std::distance(v5, v) == 1);

		auto v6 = ++v2;
		assert(v6 == v2);
		assert(v6 == v);
		assert(v2 == v);
		assert((v2 - v) == 0);
		assert(std::distance(v, v2) == 0);
		assert(std::distance(v2, v) == 0);
	}
	assert(m.end() == m.end());
	assert((m.end() + 0) == m.end());
	assert((m.end() - 0) == m.end());
	assert((m.end() - 1) != m.end());
	assert((m.end() - 1) == (m.begin() + (m.rows() - 1)));

	assert((std::size_t)std::distance(m.begin(), m.end()) == m.rows());
}
void access_tests()
{
	matrix<int> m(4, 5);

	std::size_t val = 0;
	for (std::size_t i = 0; i < m.rows(); ++i)
		for (std::size_t j = 0; j < m.cols(); ++j)
			m.at(i, j) = (int)val++;

	val = 0;
	for (std::size_t i = 0; i < m.rows(); ++i)
		for (std::size_t j = 0; j < m.cols(); ++j)
		{
			assert(m(i, j) == val++);
		}

	val = 0;
	for (auto row : m)
		for (auto v : row)
		{
			assert(v == val++);
		}

	val = 0;
	for (auto v : m.flat())
	{
		assert(v == val++);
	}

	for (auto row : m)
	{
		assert(row.size() == m.cols());
	}
	for (std::size_t i = 0; i < m.rows(); ++i)
	{
		assert(m[i].size() == m.cols());
		assert(m.at(i).size() == m.cols());
	}
	assert(m.flat().size() == m.size());
	assert(m.size() == m.rows() * m.cols());
}
void resize_flat_tests()
{
	matrix<int> m(4, 5);
	assert(m.rows() == 4);
	assert(m.cols() == 5);
	assert(m.size() == 20);

	matrix<int>::flat_view f = m.flat();
	assert(f.size() == 20);
	for (std::size_t i = 0; i < 20; ++i) f.at(i) = (int)(i + 1);
	assert_throws(f.at(20), std::out_of_range);

	m.resize_flat(5, 6, 7567);
	assert(m.rows() == 5);
	assert(m.cols() == 6);
	assert(m.size() == 30);

	f = m.flat();
	assert(f.size() == 30);
	for (std::size_t i = 0; i < 20; ++i) { assert(f.at(i) == i + 1); }
	for (std::size_t i = 20; i < 30; ++i) { assert(f.at(i) == 7567); }
	assert_throws(f.at(30), std::out_of_range);

	m.resize_flat(3, 2, 7567);
	assert(m.rows() == 3);
	assert(m.cols() == 2);
	assert(m.size() == 6);

	f = m.flat();
	assert(f.size() == 6);
	for (std::size_t i = 0; i < 6; ++i) { assert(f.at(i) == i + 1); }
	assert_throws(f.at(6), std::out_of_range);

	m.resize_flat(4, 5, 7567);
	assert(m.rows() == 4);
	assert(m.cols() == 5);
	assert(m.size() == 20);

	f = m.flat();
	assert(f.size() == 20);
	for (std::size_t i = 0; i < 6; ++i) { assert(f.at(i) == i + 1); }
	for (std::size_t i = 6; i < 20; ++i) { assert(f.at(i) == 7567); }
	assert_throws(f.at(20), std::out_of_range);

	auto cpy = m;
	assert(cpy.rows() == 4);
	assert(cpy.cols() == 5);
	assert(cpy.size() == 20);
	assert(!cpy.empty());
	cpy.resize_flat(0, 3);
	assert(cpy.rows() == 0);
	assert(cpy.cols() == 0);
	assert(cpy.size() == 0);
	assert(cpy.empty());

	cpy = m;
	assert(cpy.rows() == 4);
	assert(cpy.cols() == 5);
	assert(cpy.size() == 20);
	assert(!cpy.empty());
	cpy.resize_flat(2, 0);
	assert(cpy.rows() == 0);
	assert(cpy.cols() == 0);
	assert(cpy.size() == 0);
	assert(cpy.empty());

	cpy = m;
	assert(cpy.rows() == 4);
	assert(cpy.cols() == 5);
	assert(cpy.size() == 20);
	assert(!cpy.empty());
	cpy.resize_flat(0, 0);
	assert(cpy.rows() == 0);
	assert(cpy.cols() == 0);
	assert(cpy.size() == 0);
	assert(cpy.empty());

	cpy = m;
	assert(cpy.rows() == 4);
	assert(cpy.cols() == 5);
	assert(cpy.size() == 20);
	assert(!cpy.empty());
	assert(cpy == m);
	cpy.resize_flat(4, 5);
	assert(cpy.rows() == 4);
	assert(cpy.cols() == 5);
	assert(cpy.size() == 20);
	assert(!cpy.empty());
	assert(cpy == m);
}
void resize_cols_test()
{
	matrix<int> m(4, 5);
	assert(m.rows() == 4);
	assert(m.cols() == 5);
	assert(m.size() == 20);

	for (std::size_t i = 0; i < 4; ++i)
		for (std::size_t j = 0; j < 5; ++j)
		{
			m.at(i, j) = (int)(i * 5 + j);
		}

	{
		std::size_t val = 0;
		for (auto v : m.flat()) { assert(v == val++); }
	}

	m.resize_cols(7, 39582);
	assert(m.rows() == 4);
	assert(m.cols() == 7);
	assert(m.size() == 28);
	for (std::size_t i = 0; i < 4; ++i)
	{
		for (std::size_t j = 0; j < 5; ++j) { assert(m.at(i, j) == i * 5 + j); }
		for (std::size_t j = 5; j < 7; ++j) { assert(m.at(i, j) == 39582); }

		for (std::size_t j = 0; j < 5; ++j) { assert(m.at(i).at(j) == i * 5 + j); }
		for (std::size_t j = 5; j < 7; ++j) { assert(m.at(i).at(j) == 39582); }

		assert_throws(m.at(i, 7), std::out_of_range);
	}
	assert_throws(m.at(4), std::out_of_range);

	m.resize_cols(2, 29083495);
	assert(m.rows() == 4);
	assert(m.cols() == 2);
	assert(m.size() == 8);
	for (std::size_t i = 0; i < 4; ++i)
	{
		for (std::size_t j = 0; j < 2; ++j) { assert(m.at(i, j) == i * 5 + j); }
		assert_throws(m.at(i, 2), std::out_of_range);
	}
	assert_throws(m.at(4), std::out_of_range);

	m.resize_cols(33);
	assert(m.rows() == 4);
	assert(m.cols() == 33);
	assert(m.size() == 33 * 4);
	for (std::size_t i = 0; i < 4; ++i)
	{
		for (std::size_t j = 0; j < 2; ++j) { assert(m.at(i, j) == i * 5 + j); }
		for (std::size_t j = 2; j < 33; ++j) { assert(m.at(i, j) == 0); }

		for (std::size_t j = 0; j < 2; ++j) { assert(m.at(i).at(j) == i * 5 + j); }
		for (std::size_t j = 2; j < 33; ++j) { assert(m.at(i).at(j) == 0); }

		assert_throws(m.at(i, 33), std::out_of_range);
	}
	assert_throws(m.at(4), std::out_of_range);

	assert(!m.empty());
	m.resize_cols(0);
	assert(m.rows() == 0);
	assert(m.cols() == 0);
	assert(m.size() == 0);
	assert(m.flat().size() == 0);
	assert(m.empty());
}
void resize_rows_test()
{
	matrix<int> m(4, 5);
	assert(m.rows() == 4);
	assert(m.cols() == 5);
	assert(m.size() == 20);

	for (std::size_t i = 0; i < 4; ++i)
		for (std::size_t j = 0; j < 5; ++j)
			m.at(i, j) = (int)(i * 5 + j);

	for (std::size_t i = 0; i < 4; ++i)
	{
		for (std::size_t j = 0; j < 5; ++j)
		{
			assert(m.at(i, j) == i * 5 + j);
		}
		assert_throws(m.at(i, 5), std::out_of_range);
	}
	assert_throws(m.at(4), std::out_of_range);

	m.resize_rows(7, 204939);
	assert(m.rows() == 7);
	assert(m.cols() == 5);
	assert(m.size() == 35);

	for (std::size_t i = 0; i < 4; ++i)
	{
		for (std::size_t j = 0; j < 5; ++j)
		{
			assert(m.at(i, j) == i * 5 + j);
		}
		assert_throws(m.at(i, 5), std::out_of_range);
	}
	for (std::size_t i = 4; i < 7; ++i)
	{
		for (std::size_t j = 0; j < 5; ++j)
		{
			assert(m.at(i, j) == 204939);
		}
		assert_throws(m.at(i, 5), std::out_of_range);
	}
	assert_throws(m.at(7), std::out_of_range);

	m.resize_rows(2);
	assert(m.rows() == 2);
	assert(m.cols() == 5);
	assert(m.size() == 10);

	for (std::size_t i = 0; i < 2; ++i)
	{
		for (std::size_t j = 0; j < 5; ++j)
		{
			assert(m.at(i, j) == i * 5 + j);
		}
		assert_throws(m.at(i, 5), std::out_of_range);
	}
	assert_throws(m.at(2), std::out_of_range);

	m.resize_rows(4);
	assert(m.rows() == 4);
	assert(m.cols() == 5);
	assert(m.size() == 20);

	for (std::size_t i = 0; i < 2; ++i)
	{
		for (std::size_t j = 0; j < 5; ++j)
		{
			assert(m.at(i, j) == i * 5 + j);
		}
		assert_throws(m.at(i, 5), std::out_of_range);
	}
	for (std::size_t i = 2; i < 4; ++i)
	{
		for (std::size_t j = 0; j < 5; ++j)
		{
			assert(m.at(i, j) == 0);
		}
		assert_throws(m.at(i, 5), std::out_of_range);
	}
	assert_throws(m.at(4), std::out_of_range);

	assert(!m.empty());
	m.resize_rows(0);
	assert(m.rows() == 0);
	assert(m.cols() == 0);
	assert(m.size() == 0);
	assert(m.empty());
}
void resize_tests()
{
	matrix<int> m(4, 5);
	assert(m.rows() == 4);
	assert(m.cols() == 5);
	assert(m.size() == 20);
	assert(!m.empty());

	for (std::size_t i = 0; i < 4; ++i)
		for (std::size_t j = 0; j < 5; ++j)
			m(i, j) = (int)(i * 5 + j);

	m.resize(3, 8, 3574835);
	assert(m.rows() == 3);
	assert(m.cols() == 8);
	assert(m.size() == 24);
	assert(!m.empty());

	for (std::size_t i = 0; i < 3; ++i)
	{
		for (std::size_t j = 0; j < 5; ++j) { assert(m.at(i, j) == i * 5 + j); }
		for (std::size_t j = 5; j < 8; ++j) { assert(m.at(i, j) == 3574835); }
		assert_throws(m.at(i, 8), std::out_of_range);
	}
	assert_throws(m.at(3), std::out_of_range);

	m.resize(7, 3, -57754);
	assert(m.rows() == 7);
	assert(m.cols() == 3);
	assert(m.size() == 21);
	assert(!m.empty());

	for (std::size_t i = 0; i < 3; ++i)
	{
		for (std::size_t j = 0; j < 3; ++j) { assert(m.at(i, j) == i * 5 + j); }
		assert_throws(m.at(i, 3), std::out_of_range);
	}
	for (std::size_t i = 3; i < 7; ++i)
	{
		for (std::size_t j = 0; j < 3; ++j) { assert(m.at(i, j) == -57754); }
		assert_throws(m.at(i, 3), std::out_of_range);
	}
	assert_throws(m.at(7), std::out_of_range);

	m.resize(9, 11);
	assert(m.rows() == 9);
	assert(m.cols() == 11);
	assert(m.size() == 99);
	assert(!m.empty());

	for (std::size_t i = 0; i < 3; ++i)
	{
		for (std::size_t j = 0; j < 3; ++j) { assert(m.at(i, j) == i * 5 + j); }
		for (std::size_t j = 3; j < 11; ++j) { assert(m.at(i, j) == 0); }
		assert_throws(m.at(i, 11), std::out_of_range);
	}
	for (std::size_t i = 3; i < 7; ++i)
	{
		for (std::size_t j = 0; j < 3; ++j) { assert(m.at(i, j) == -57754); }
		for (std::size_t j = 3; j < 11; ++j) { assert(m.at(i, j) == 0); }
		assert_throws(m.at(i, 11), std::out_of_range);
	}
	for (std::size_t i = 7; i < 9; ++i)
	{
		for (std::size_t j = 0; j < 11; ++j) { assert(m.at(i, j) == 0); }
		assert_throws(m.at(i, 11), std::out_of_range);
	}
	assert_throws(m.at(9), std::out_of_range);

	m.resize(5, 8);
	assert(m.rows() == 5);
	assert(m.cols() == 8);
	assert(m.size() == 40);
	assert(!m.empty());

	for (std::size_t i = 0; i < 3; ++i)
	{
		for (std::size_t j = 0; j < 3; ++j) { assert(m.at(i, j) == i * 5 + j); }
		for (std::size_t j = 3; j < 8; ++j) { assert(m.at(i, j) == 0); }
		assert_throws(m.at(i, 8), std::out_of_range);
	}
	for (std::size_t i = 3; i < 5; ++i)
	{
		for (std::size_t j = 0; j < 3; ++j) { assert(m.at(i, j) == -57754); }
		for (std::size_t j = 3; j < 8; ++j) { assert(m.at(i, j) == 0); }
		assert_throws(m.at(i, 8), std::out_of_range);
	}
	assert_throws(m.at(5), std::out_of_range);

	auto cpy = m;
	assert(m.rows() == 5);
	assert(m.cols() == 8);
	assert(m.size() == 40);
	assert(!m.empty());
	assert(cpy.rows() == 5);
	assert(cpy.cols() == 8);
	assert(cpy.size() == 40);
	assert(!cpy.empty());
	assert(cpy == m);

	cpy.resize(0, 43);
	assert(cpy.rows() == 0);
	assert(cpy.cols() == 0);
	assert(cpy.size() == 0);
	assert(cpy.flat().size() == 0);
	assert(cpy.empty());
	assert(cpy != m);

	cpy = m;
	assert(m.rows() == 5);
	assert(m.cols() == 8);
	assert(m.size() == 40);
	assert(!m.empty());
	assert(cpy.rows() == 5);
	assert(cpy.cols() == 8);
	assert(cpy.size() == 40);
	assert(!cpy.empty());
	assert(cpy == m);

	cpy.resize(14, 0);
	assert(cpy.rows() == 0);
	assert(cpy.cols() == 0);
	assert(cpy.size() == 0);
	assert(cpy.flat().size() == 0);
	assert(cpy.empty());
	assert(cpy != m);
}
void append_cols_tests()
{
	typedef intrinsic<int> t;

	matrix<t> m = { { (t)1, (t)2, (t)3 }, { (t)4, (t)5, (t)6 }, { (t)7, (t)8, (t)9 }, { (t)10, (t)11, (t)12 } };
	matrix<t> n = { { (t)77, (t)64 }, { (t)100, (t)(-123) }, { (t)54, (t)65 }, { (t)7, (t)0 } };

	assert(m.rows() == 4);
	assert(m.cols() == 3);
	assert(m.size() == 12);
	assert(!m.empty());

	assert(n.rows() == 4);
	assert(n.cols() == 2);
	assert(n.size() == 8);
	assert(!n.empty());

	m.append_cols(n);
	
	assert(m.rows() == 4);
	assert(m.cols() == 5);
	assert(m.size() == 20);
	assert(!m.empty());

	assert(n.rows() == 4);
	assert(n.cols() == 2);
	assert(n.size() == 8);
	assert(!n.empty());

	assert((m == matrix<t>{ { (t)1, (t)2, (t)3, (t)77, (t)64 }, { (t)4, (t)5, (t)6, (t)100, (t)(-123) }, { (t)7, (t)8, (t)9, (t)54, (t)65 }, { (t)10, (t)11, (t)12, (t)7, (t)0 } }));
	assert((n == matrix<t>{ { (t)77, (t)64 }, { (t)100, (t)(-123) }, { (t)54, (t)65 }, { (t)7, (t)0 } }));

	n.append_cols(std::move(m));

	assert(m.rows() == 0);
	assert(m.cols() == 0);
	assert(m.size() == 0);
	assert(m.empty());

	assert(n.rows() == 4);
	assert(n.cols() == 7);
	assert(n.size() == 28);
	assert(!n.empty());

	assert((m == matrix<t>{ }));
	assert((n == matrix<t>{ { (t)77, (t)64, (t)1, (t)2, (t)3, (t)77, (t)64 }, { (t)100, (t)(-123), (t)4, (t)5, (t)6, (t)100, (t)(-123) },
		{ (t)54, (t)65, (t)7, (t)8, (t)9, (t)54, (t)65 }, { (t)7, (t)0, (t)10, (t)11, (t)12, (t)7, (t)0 } }));

	assert_throws(n.append_cols(m), matrix_size_error);
	assert_throws(m.append_cols(n), matrix_size_error);

	n.append_cols(n);
	assert(n.rows() == 4);
	assert(n.cols() == 14);
	assert(n.size() == 56);
	assert(!n.empty());

	assert((n == matrix<t>{ { (t)77, (t)64, (t)1, (t)2, (t)3, (t)77, (t)64, (t)77, (t)64, (t)1, (t)2, (t)3, (t)77, (t)64 },
		{ (t)100, (t)(-123), (t)4, (t)5, (t)6, (t)100, (t)(-123), (t)100, (t)(-123), (t)4, (t)5, (t)6, (t)100, (t)(-123) },
		{ (t)54, (t)65, (t)7, (t)8, (t)9, (t)54, (t)65, (t)54, (t)65, (t)7, (t)8, (t)9, (t)54, (t)65 },
		{ (t)7, (t)0, (t)10, (t)11, (t)12, (t)7, (t)0, (t)7, (t)0, (t)10, (t)11, (t)12, (t)7, (t)0 } }));
}

void comparison_tests()
{
	matrix<int> m(4, 5);
	matrix<int> n = m;

	assert(m == n);
	assert(!(m != n));

	m.resize_flat(5, 4);
	{
		const auto f1 = m.flat();
		const auto f2 = n.flat();
		assert(f1.size() == f2.size());
		for (std::size_t i = 0; i < f1.size(); ++i) { assert(f1[i] == f2[i]); }
	}
	assert(m != n);
	assert(!(m == n));

	m = n;
	assert(m == n);
	assert(!(m != n));

	m.at(2, 2) = 5;
	assert(m != n);
	assert(!(m == n));

	n.at(2, 2) = 5;
	assert(m == n);
	assert(!(m != n));
}

void add_tests()
{
	matrix<std::size_t> a(4, 5), b, res;
	assert(b.empty());
	assert(res.empty());

	b.resize_flat(a.rows(), a.cols() - 1);
	assert_throws(a + b, matrix_size_error);
	b.resize_flat(a.rows(), a.cols() + 1);
	assert_throws(a + b, matrix_size_error);

	b.resize_flat(a.rows() - 1, a.cols());
	assert_throws(a + b, matrix_size_error);
	b.resize_flat(a.rows() + 1, a.cols());
	assert_throws(a + b, matrix_size_error);

	b.resize_flat(a.rows(), a.cols());
	for (std::size_t i = 0; i < a.rows(); ++i)
		for (std::size_t j = 0; j < a.cols(); ++j)
		{
			a.at(i, j) = i + j * j + 18;
			b.at(i, j) = i * i + 12;
		}
	assert(a != b);

	res = a + b;
	assert(!res.empty());
	assert(res.rows() == a.rows());
	assert(res.cols() == a.cols());

	for (std::size_t i = 0; i < res.rows(); ++i)
		for (std::size_t j = 0; j < res.cols(); ++j)
		{
			assert(res.at(i, j) == i + j * j + 18 + i * i + 12);
		}
}
void sub_tests()
{
	matrix<int> a(5, 7), b, res;

	b.resize_flat(a.rows(), a.cols() - 1);
	assert_throws(a - b, matrix_size_error);
	b.resize_flat(a.rows(), a.cols() + 1);
	assert_throws(a - b, matrix_size_error);

	b.resize_flat(a.rows() - 1, a.cols());
	assert_throws(a - b, matrix_size_error);
	b.resize_flat(a.rows() + 1, a.cols());
	assert_throws(a - b, matrix_size_error);

	b.resize_flat(a.rows(), a.cols());
	for (std::size_t i = 0; i < a.rows(); ++i)
		for (std::size_t j = 0; j < a.cols(); ++j)
		{
			a.at(i, j) = (int)to_signed(2 * i + j * j);
			b(i, j) = (int)to_signed(j + j * j - i);
		}

	res = a - b;
	assert(!res.empty());
	assert(res.rows() == a.rows());
	assert(res.cols() == a.cols());

	for (std::size_t i = 0; i < a.rows(); ++i)
		for (std::size_t j = 0; j < a.cols(); ++j)
		{
			assert(res.at(i, j) == (int)to_signed((2 * i + j * j) - (j + j * j - i)));
		}
}
void mul_tests()
{
	matrix<int> a = { { 1, 2, 3 }, { 4, 5, 6 } };
	matrix<int> b = { { 7, 8 }, { 9, 10 }, { 11, 12 } };
	matrix<int> c = a * b;
	assert((c == matrix<int>{ { 58, 64 }, { 139, 154 } }));
	c = b * a;
	assert((c == matrix<int>{ { 39, 54, 69 }, { 49, 68, 87 }, { 59, 82, 105 } }));

	a = { { 3, 4, 2 } };
	b = { { 13, 9, 7, 15 }, { 8, 7, 4, 6 }, { 6, 4, 0, 3 } };
	c = a * b;
	assert((c == matrix<int>{ { 83, 63, 37, 75 } }));
	assert_throws(b * a, matrix_size_error);

	a = { { 1, 2, 3 } };
	b = { { 4 }, { 5 }, { 6 } };
	c = a * b;
	assert((c == matrix<int>{ { 32 } }));
	c = b * a;
	assert((c == matrix<int>{ { 4, 8, 12 }, { 5, 10, 15 }, { 6, 12, 18 } }));

	a = { { 1, 2 }, { 3, 4 } };
	b = { { 2, 0 }, { 1, 2 } };
	c = a * b;
	assert((c == matrix<int>{ { 4, 4 }, { 10, 8 } }));
	c = b * a;
	assert((c == matrix<int>{ { 2, 4 }, { 7, 10 } }));

	a = {};
	b.clear();
	assert(a == b);
	c = a * b;
	assert(c == a && c.empty());
	assert(b * a == c);
}

void unseq_tests()
{
	static_assert(matrix_impl::is_unseq<short>::value, "primitive unseq error");
	static_assert(matrix_impl::is_unseq<int>::value, "primitive unseq error");
	static_assert(matrix_impl::is_unseq<long>::value, "primitive unseq error");
	static_assert(matrix_impl::is_unseq<long long>::value, "primitive unseq error");

	static_assert(matrix_impl::is_unseq<unsigned short>::value, "primitive unseq error");
	static_assert(matrix_impl::is_unseq<unsigned int>::value, "primitive unseq error");
	static_assert(matrix_impl::is_unseq<unsigned long>::value, "primitive unseq error");
	static_assert(matrix_impl::is_unseq<unsigned long long>::value, "primitive unseq error");

	static_assert(matrix_impl::is_unseq<float>::value, "primitive unseq error");
	static_assert(matrix_impl::is_unseq<double>::value, "primitive unseq error");
	static_assert(matrix_impl::is_unseq<long double>::value, "primitive unseq error");

	static_assert(matrix_impl::is_unseq<matrix<int>::row_view>::value, "expr unseq error");
	static_assert(matrix_impl::is_unseq<const matrix<int>::row_view>::value, "expr unseq error");
	static_assert(matrix_impl::is_unseq<matrix<int>::row_view&>::value, "expr unseq error");
	static_assert(matrix_impl::is_unseq<const matrix<int>::row_view&>::value, "expr unseq error");
	static_assert(matrix_impl::is_unseq<matrix<int>::row_view&&>::value, "expr unseq error");
	static_assert(matrix_impl::is_unseq<const matrix<int>::row_view&&>::value, "expr unseq error");

	static_assert(!matrix_impl::is_unseq<intrinsic<int>>::value, "expr unseq error");
	static_assert(!matrix_impl::is_unseq<const intrinsic<int>>::value, "expr unseq error");
	static_assert(!matrix_impl::is_unseq<intrinsic<int>&>::value, "expr unseq error");
	static_assert(!matrix_impl::is_unseq<const intrinsic<int>&>::value, "expr unseq error");
	static_assert(!matrix_impl::is_unseq<intrinsic<int>&&>::value, "expr unseq error");
	static_assert(!matrix_impl::is_unseq<const intrinsic<int>&&>::value, "expr unseq error");

	{
		matrix<int> m;
		static_assert(matrix_impl::is_unseq<decltype(m[0])>::value, "expr unseq error");
		static_assert(matrix_impl::is_unseq<decltype(+m[0])>::value, "expr unseq error");
		static_assert(matrix_impl::is_unseq<decltype(-m[0])>::value, "expr unseq error");
		static_assert(matrix_impl::is_unseq<decltype(m[0] + m[0])>::value, "expr unseq error");
		static_assert(matrix_impl::is_unseq<decltype(m[0] - m[0])>::value, "expr unseq error");
		static_assert(matrix_impl::is_unseq<decltype(m[0] - -m[0])>::value, "expr unseq error");
		static_assert(matrix_impl::is_unseq<decltype(m[0] - -+m[0])>::value, "expr unseq error");
		static_assert(matrix_impl::is_unseq<decltype(+-m[0] - -+m[0])>::value, "expr unseq error");
		static_assert(matrix_impl::is_unseq<decltype(+-m[0] * 3 - -+m[0])>::value, "expr unseq error");
		static_assert(matrix_impl::is_unseq<decltype(3 * +-m[0] - -+m[0])>::value, "expr unseq error");
	}
	{
		matrix<double> m;
		static_assert(matrix_impl::is_unseq<decltype(m[0])>::value, "expr unseq error");
		static_assert(matrix_impl::is_unseq<decltype(+m[0])>::value, "expr unseq error");
		static_assert(matrix_impl::is_unseq<decltype(-m[0])>::value, "expr unseq error");
		static_assert(matrix_impl::is_unseq<decltype(m[0] + m[0])>::value, "expr unseq error");
		static_assert(matrix_impl::is_unseq<decltype(m[0] - m[0])>::value, "expr unseq error");
		static_assert(matrix_impl::is_unseq<decltype(m[0] - -m[0])>::value, "expr unseq error");
		static_assert(matrix_impl::is_unseq<decltype(m[0] - -+m[0])>::value, "expr unseq error");
		static_assert(matrix_impl::is_unseq<decltype(+-m[0] - -+m[0])>::value, "expr unseq error");
		static_assert(matrix_impl::is_unseq<decltype(+-m[0] * 3 - -+m[0])>::value, "expr unseq error");
		static_assert(matrix_impl::is_unseq<decltype(3 * +-m[0] - -+m[0])>::value, "expr unseq error");
	}
	{
		matrix<intrinsic<int>> m;
		static_assert(!matrix_impl::is_unseq<decltype(m[0])>::value, "expr unseq error");
		static_assert(!matrix_impl::is_unseq<decltype(+m[0])>::value, "expr unseq error");
		static_assert(!matrix_impl::is_unseq<decltype(-m[0])>::value, "expr unseq error");
		static_assert(!matrix_impl::is_unseq<decltype(m[0] * 3)>::value, "expr unseq error");
		static_assert(!matrix_impl::is_unseq<decltype(3 * m[0])>::value, "expr unseq error");
		static_assert(!matrix_impl::is_unseq<decltype(m[0] + m[0])>::value, "expr unseq error");
		static_assert(!matrix_impl::is_unseq<decltype(m[0] - m[0])>::value, "expr unseq error");
		static_assert(!matrix_impl::is_unseq<decltype(m[0] - -m[0])>::value, "expr unseq error");
		static_assert(!matrix_impl::is_unseq<decltype(m[0] - -+m[0])>::value, "expr unseq error");
		static_assert(!matrix_impl::is_unseq<decltype(+-m[0] - -+m[0])>::value, "expr unseq error");
		static_assert(!matrix_impl::is_unseq<decltype(+-m[0] * 3 - -+m[0])>::value, "expr unseq error");
		static_assert(!matrix_impl::is_unseq<decltype(3 * +-m[0] - -+m[0])>::value, "expr unseq error");
		static_assert(!matrix_impl::is_unseq<decltype(3 * -m[0])>::value, "expr unseq error");
		static_assert(!matrix_impl::is_unseq<decltype(intrinsic<int>{3} * -m[0]) > ::value, "expr unseq error");
	}
	{
		typedef intrinsic<int> t;
		static_assert(std::is_assignable_v<int&, t>, "test assignability failure");

		matrix<int> m1 = { { 5, 2 }, { 4, 1 }, { 6, 7 } };
		matrix<intrinsic<int>> m2 = { { (t)1, (t)2 }, { (t)3, (t)4 } };

		m1[0] = m1[0] + m2[1];
		assert((m1 == matrix<int>{ { 8, 6 }, { 4, 1 }, { 6, 7 } }));

		static_assert(matrix_impl::is_unseq<decltype(m1[0])>::value, "expr unseq error");
		static_assert(!matrix_impl::is_unseq<decltype(m2[0])>::value, "expr unseq error");
		static_assert(!matrix_impl::is_unseq<decltype(m1[0] + m2[0])>::value, "expr unseq error");
		static_assert(!matrix_impl::is_unseq<decltype(m1[0] + m2[0])>::value, "expr unseq error");
	}
}
void row_op_tests()
{
	static_assert(!std::is_assignable_v<matrix<int>::row_view&&, matrix<int>::row_view>, "this is disallowed to avoid being misleading");

	matrix<int> v = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };

	v[0] = v[1] + v[2];
	assert((v == matrix<int>{ {11, 13, 15}, { 4, 5, 6 }, { 7, 8, 9 } }));

	v[1] = +v[0]; // without the + this would fail to compile (by design - guaranteed by above assertion)
	assert((v == matrix<int>{ {11, 13, 15}, { 11, 13, 15 }, { 7, 8, 9 } }));

	v[0] = v[1] - v[2];
	assert((v == matrix<int>{ {4, 5, 6}, { 11, 13, 15 }, { 7, 8, 9 } }));

	v[2] = -v[1];
	assert((v == matrix<int>{ {4, 5, 6}, { 11, 13, 15 }, { -11, -13, -15 } }));

	v[1] = v[0] * 2;
	assert((v == matrix<int>{ {4, 5, 6}, { 8, 10, 12 }, { -11, -13, -15 } }));

	v[0] = -3 * v[2];
	assert((v == matrix<int>{ {33, 39, 45}, { 8, 10, 12 }, { -11, -13, -15 } }));

	v[1] = v[1] / 2;
	assert((v == matrix<int>{ {33, 39, 45}, { 4, 5, 6 }, { -11, -13, -15 } }));

	v[0][0] = 3;
	assert((v == matrix<int>{ {3, 39, 45}, { 4, 5, 6 }, { -11, -13, -15 } }));

	v[0] = v[0] / v[0][0];
	assert((v == matrix<int>{ {1, 13, 15}, { 4, 5, 6 }, { -11, -13, -15 } }));

	{
		static_assert(std::is_same_v<decltype(v[0][0]), int&>, "wrong type returned");
		static_assert(std::is_same_v<decltype(v[0].at(0)), int&>, "wrong type returned");
		static_assert(std::is_same_v<decltype(v[0].move()[0]), int&&>, "wrong type returned");

		const auto &m = v;
		static_assert(std::is_same_v<decltype(m[0][0]), const int&>, "wrong type returned");
		static_assert(std::is_same_v<decltype(m[0].at(0)), const int&>, "wrong type returned");
	}

	v[1] = 3 * v[0] - v[2] / v[1][1];
	assert((v == matrix<int>{ {1, 13, 15}, { 5, 41, 48 }, { -11, -13, -15 } }));

	v[0] = 2 * -v[2];
	assert((v == matrix<int>{ {22, 26, 30}, { 5, 41, 48 }, { -11, -13, -15 } }));

	v[0] -= v[1];
	assert((v == matrix<int>{ {17, -15, -18}, { 5, 41, 48 }, { -11, -13, -15 } }));

	v[2] -= 2 * v[0];
	assert((v == matrix<int>{ {17, -15, -18}, { 5, 41, 48 }, { -45, 17, 21 } }));

	v[1] += v[0];
	assert((v == matrix<int>{ {17, -15, -18}, { 22, 26, 30 }, { -45, 17, 21 } }));

	v[0] += 3 * v[2];
	assert((v == matrix<int>{ {-118, 36, 45}, { 22, 26, 30 }, { -45, 17, 21 } }));

	v[0] += v[0];
	assert((v == matrix<int>{ {-236, 72, 90}, { 22, 26, 30 }, { -45, 17, 21 } }));

	v[0] -= v[0];
	assert((v == matrix<int>{ {0, 0, 0}, { 22, 26, 30 }, { -45, 17, 21 } }));

	{
		typedef intrinsic<int> t;
		matrix<t> m{ { (t)1, (t)3, (t)5 }, { (t)2, (t)4, (t)3 }, { (t)6, (t)2, (t)4 } };

		m[1] = +m[0];
		assert((m == matrix<t>{ { (t)1, (t)3, (t)5 }, { (t)1, (t)3, (t)5 }, { (t)6, (t)2, (t)4 } }));

		m[0] = m[2].move();
		assert((m == matrix<t>{ { (t)6, (t)2, (t)4 }, { (t)1, (t)3, (t)5 }, { (t)t::move_assign, (t)t::move_assign, (t)t::move_assign } }));
	}
}

int main() try
{
	std::cout << "beginning tests\n" << std::flush;

	meta_tests();

	basic_tests();
	init_list_tests();
	iter_tests();
	access_tests();
	resize_flat_tests();
	resize_cols_test();
	resize_rows_test();
	resize_tests();
	append_cols_tests();

	comparison_tests();

	add_tests();
	sub_tests();
	mul_tests();

	unseq_tests();
	row_op_tests();

	std::cout << "all tests passed successfully\n\nbeginning benchmarks\n" << std::flush;

	benchmark_binary<double>("copy", 100, { 1000, 1000 }, { 0, 0 }, [](const auto &a, const auto &b) { return a; });

	benchmark_binary<double>("add assign", 100, { 1000, 1000 }, { 1000, 1000 }, [](auto &a, const auto &b) -> auto& { return a += b; });
	benchmark_binary<double>("add", 100, { 1000, 1000 }, { 1000, 1000 }, [](const auto &a, const auto &b) { return a + b; });

	benchmark_binary<double>("sub assign", 100, { 1000, 1000 }, { 1000, 1000 }, [](auto &a, const auto &b) -> auto& { return a -= b; });
	benchmark_binary<double>("sub", 100, { 1000, 1000 }, { 1000, 1000 }, [](const auto &a, const auto &b) { return a - b; });

	benchmark_binary<double>("mul matrix", 10, { 1000, 1000 }, { 1000, 1000 }, [](const auto &a, const auto &b) { return a * b; });

	benchmark_binary<double>("mul scalar assign", 100, { 1000, 1000 }, { 1, 1 }, [](auto &a, const auto &b) -> auto& { return a *= b.flat()[0]; });
	benchmark_binary<double>("mul scalar", 100, { 1000, 1000 }, { 1, 1 }, [](const auto &a, const auto &b) { return a * b.flat()[0]; });

	std::cout << "all benchmarks completed\n" << std::flush;
	return 0;
}
catch (const std::exception & ex)
{
	std::cerr << "\n\nUNHANDLED EXCEPTION: " << ex.what() << '\n';
}
catch (...)
{
	std::cerr << "\n\nUNHANDLED EXCEPTION OF UNKNOWN TYPE\n";
}
