#include <iostream>
#include <iomanip>
#include "Matricies.h"

#define forall(m) for(std::size_t row = 0; row < m.rows(); ++row) for(std::size_t col = 0; col < m.cols(); ++col)

int main()
{
	try
	{
		Matrix<double> m;
		std::size_t rows, cols;
		double det;

		while (true)
		{
			std::cout << "\n\nrows: "; std::cin >> rows;
			std::cout << "cols: "; std::cin >> cols;

			if (std::cin.fail()) { std::cin.clear(); std::cin.ignore(32767, '\n'); continue; }

			m.resize_dump(rows, cols);
			std::cout << "data: ";
			for (std::size_t i = 0; i < rows * cols; ++i) std::cin >> m[i];

			if (std::cin.fail()) { std::cin.clear(); std::cin.ignore(32767, '\n'); continue; }

			std::cout << '\n' << m;
			std::cout << "rank: " << m.RREF(&det) << '\n';
			std::cout << "det:  " << det << '\n';
			std::cout << m << '\n';
		}
	}
	catch (const std::exception &ex)
	{
		std::cerr << "error: " << ex.what() << '\n';
	}

	std::cin.get();
	return 0;
}