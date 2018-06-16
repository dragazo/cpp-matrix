#include <iostream>
#include <iomanip>
#include <chrono>
#include "Matricies.h"

#define forall(m) for(std::size_t row = 0; row < m.rows(); ++row) for(std::size_t col = 0; col < m.cols(); ++col)

#define __BENCHMARK 0

int main()
{
	using namespace std::chrono;
	typedef high_resolution_clock hrc;

	try
	{
		#if __BENCHMARK
		Matrix<float> a(10, 10), b(10, 10);

		forall(a)
		{
			a(row, col) = row + col + row * col;
			b(row, col) = row * row + col * col;
		}

		const int reps = 20;
		
		auto start = hrc::now();

		for (int i = 0; i < reps; ++i)
		{
			a *= 1.1f;
		}

		auto t = duration_cast<milliseconds>(hrc::now() - start).count();

		//std::cout << a << "\n\n\n";

		std::cout << "dim:  " << a.rows() << " x " << a.cols() << '\n';
		std::cout << "reps: " << reps << '\n';
		std::cout << "elapsed: " << t << "ms\n";
		std::cout << "average: " << std::setprecision(3) << ((double)t / reps) << "ms\n";

		std::cout << "\n\n\n" << a << '\n';
		std::cin.get();

		#else

		Matrix<double> m, temp;
		std::size_t rows, cols;
		double det;
		
		while (true)
		{
			std::cout << "rows: "; std::cin >> rows;
			std::cout << "cols: "; std::cin >> cols;

			if (std::cin.fail()) { std::cin.clear(); std::cin.ignore(32767, '\n'); continue; }

			m.resize_dump(rows, cols);
			std::cout << "data: ";
			for (std::size_t i = 0; i < rows * cols; ++i) std::cin >> m[i];

			if (std::cin.fail()) { std::cin.clear(); std::cin.ignore(32767, '\n'); continue; }

			temp = m;

			std::cout << '\n' << m;
			std::cout << "rank: " << m.RREF(&det) << '\n';
			std::cout << "det:  " << det << '\n';
			std::cout << m << '\n';
			std::cout << "conj:\n" << temp.conj() << '\n';

			std::cout << "\n\n";
		}
		
		#endif
	}
	catch (const std::exception &ex)
	{
		std::cerr << "error: " << ex.what() << '\n';
	}

	std::cin.get();
	return 0;
}