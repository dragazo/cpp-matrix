#include <iostream>
#include <iomanip>
#include "Matricies.h"

int main()
{
	try
	{
		Matrix<double> m(3, 4);
		Matrix<double> other;

		m(0, 0) = 1; m(0, 1) = 5; m(0, 2) = 0; m(0, 3) = 7;
		m(1, 0) = 1; m(1, 1) = 2; m(1, 2) = 1; m(1, 3) = 7;
		m(2, 0) = 5; m(2, 1) = 0; m(2, 2) = 1; m(2, 3) = 0;
		//m(3, 0) = 1; m(3, 1) = 2; m(3, 2) = 1; m(3, 3) = 8;

		std::cout << m << '\n';

		//std::cout << "det:      " << m.det() << '\n';
		std::cout << "rank:     " << m.RREF() << '\n' << m << '\n';
	}
	catch (const std::exception &ex)
	{
		std::cerr << "error: " << ex.what() << '\n';
	}

	std::cin.get();
	return 0;
}