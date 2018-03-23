#ifndef MATRICIES_H
#define MATRICIES_H

#include <cstdlib>
#include <vector>

template<typename T>
class Matrix
{
public: // -- enums / etc -- /

	enum OperationError
	{
		SizeError // operation failed due to incompatible size(s)
	};

private: // -- data -- //

	std::vector<T> data; // the elements in the array
	std::size_t r, c;    // number of rows / cols

public: // -- ctor / dtor / asgn -- //

	// due to using std::vector for the data, default ctor / dtor / asgn are sufficient

public: // -- operations -- //

};

#endif
