#ifndef C_READER_H
#define C_READER_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>     // std::istream_iterator
#include "../DataColle/cgal_def.h"
#include <list>

class CReader
{
public:
	static void read_data_array_file(std::string filename, std::vector<Point_2> & pts, std::vector<bool> & is_constrained_list);
};


#endif // !C_READER_H


