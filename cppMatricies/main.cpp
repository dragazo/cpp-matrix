#include <iostream>
#include <iomanip>
#include "Matricies.h"



int main()
{
	Matrix<double> m(3, 4);

	m(0, 0) = 1; m(0, 1) = 5; m(0, 2) = 0; m(0, 3) = 7;
	m(1, 0) = 1; m(1, 1) = 2; m(1, 2) = 1; m(1, 3) = 7;
	m(2, 0) = 5; m(2, 1) = 0; m(2, 2) = 1; m(2, 3) = 0;

	std::cout << m << '\n';

	m.REF();
	//m.divRow(1, 2);

	std::cout << m << '\n';

	std::cin.get();
	return 0;
}