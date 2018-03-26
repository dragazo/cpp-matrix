#include <iostream>
#include <iomanip>
#include "Matricies.h"

#define forall(m) for(std::size_t row = 0; row < m.rows(); ++row) for(std::size_t col = 0; col < m.cols(); ++col)

int main()
{
	try
	{
		Matrix<double> a(3, 3);
		Matrix<double> b(3, 3);

		forall(a)
		{
			a(row, col) = row + col;
			b(row, col) = (row + col) * 20;
		}

		std::cout << a << '\n' << b << '\n';

		std::cout << "\n---------\n\n";

		a.cat_rows(b);

		std::cout << a << '\n';
	}
	catch (const std::exception &ex)
	{
		std::cerr << "error: " << ex.what() << '\n';
	}

	std::cin.get();
	return 0;
}