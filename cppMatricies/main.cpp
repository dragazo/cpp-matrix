#include <iostream>
#include <iomanip>
#include "Matricies.h"

#define forall(m) for(std::size_t row = 0; row < m.rows(); ++row) for(std::size_t col = 0; col < m.cols(); ++col)

int main()
{
	try
	{
		Matrix<double> a(100,100);
		//Matrix<double> b(3, 3);

		forall(a)
		{
			a(row, col) = std::tan(100 * row + col);
			//b(row, col) = (row + col) * 20;
		}

		std::cout << a << '\n';// << b << '\n';

		std::cout << "\n---------\n\n";

		double det;

		std::cout << "rank: " << a.REF(&det) << "\ndet: ";
		std::cout << det << '\n' << a << '\n';
	}
	catch (const std::exception &ex)
	{
		std::cerr << "error: " << ex.what() << '\n';
	}

	std::cin.get();
	return 0;
}