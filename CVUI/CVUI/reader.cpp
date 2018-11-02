#include "reader.h"
//#include "my_data_type.h"

void CReader::read_data_array_file(std::string filename, std::vector<Point_2> & pts, std::vector<bool> & is_constrained_list)
{
	pts.clear();

	std::ifstream in(filename);
	std::istream_iterator<double> iit(in);
	std::istream_iterator<double> eos;
	while (iit != eos) {
		double x, y;
		bool fixed = true;
		x = *iit; ++iit;
		y = *iit; ++iit;
		if (*iit++ == 1) fixed = false;
		
		is_constrained_list.push_back(fixed);
		pts.push_back(Point_2(x, y));
	}
}
