#include <iostream>
#include <iomanip>
#include "Matricies.h"

#define forall(m) for(std::size_t row = 0; row < m.rows(); ++row) for(std::size_t col = 0; col < m.cols(); ++col)

int main()
{
	try
	{
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

			/*
			m.resize(3, 3);
			forall(m) m(row, col) = (row + col + 1) + row * col;

			std::cout << m << '\n';

			m.resize(4, 4); std::cout << m << '\n';
			m.resize(5, 5); std::cout << m << '\n';
			m.resize(2, 2); std::cout << m << '\n';
			m.resize(6, 6); std::cout << m << '\n';

			std::cin.get();
			*/
		}
	}
	catch (const std::exception &ex)
	{
		std::cerr << "error: " << ex.what() << '\n';
	}

	std::cin.get();
	return 0;
}