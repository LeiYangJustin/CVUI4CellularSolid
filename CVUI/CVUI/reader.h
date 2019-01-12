#ifndef C_READER_H
#define C_READER_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>     // std::istream_iterator
#include "../DataColle/cgal_def.h"
#include <list>
#include "../DataColle/voronoi_diagram.h"

class CReader
{
public:
	static void read_data_array_file(std::string filename, std::vector<Point_2> & pts, std::vector<bool> & is_constrained_list);

	static bool write_triangles_from_rt(const std::string &fname, Regular_triangulation *rt);
	static bool write_vertices_nearest(const std::string &fname, Regular_triangulation *rt);
	static bool write_faces_nearest(const std::string &fname, Regular_triangulation *rt);
	static bool write_pts(const std::string &fname, std::set<Point_2> pts);
};


#endif // !C_READER_H


