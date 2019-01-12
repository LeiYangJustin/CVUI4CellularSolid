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

bool CReader::write_triangles_from_rt(const std::string & fname, Regular_triangulation * rt)
{
	std::cout << "write " << fname << std::endl;

	std::ofstream out_file;
	out_file.open(fname);
	if (out_file.is_open())
	{
		int cnt = 0;
		for (auto it = rt->finite_faces_begin();
			it != rt->finite_faces_end(); ++it)
		{
			out_file 
				<< it->vertex(0)->point().point().x()
				<< " " << it->vertex(0)->point().point().y()
				<< " " << it->vertex(1)->point().point().x()
				<< " " << it->vertex(1)->point().point().y()
				<< " " << it->vertex(2)->point().point().x()
				<< " " << it->vertex(2)->point().point().y()
				<< std::endl;
			//std::cout << ++cnt << std::endl;
		}
		out_file.close();
		return true;
	}
	return false;
}

bool CReader::write_vertices_nearest(const std::string & fname, Regular_triangulation *rt)
{
	std::cout << "write " << fname << std::endl;

	std::ofstream out_file;
	out_file.open(fname);
	if (out_file.is_open())
	{
		int cnt = 0;
		for (auto it = rt->finite_vertices_begin();
			it != rt->finite_vertices_end(); ++it)
		{
			out_file
				<< it->point().point() << " "
				<< it->nearest()
				<< std::endl;
			//std::cout << ++cnt << std::endl;
		}
		out_file.close();
		return true;
	}
	return false;
}

bool CReader::write_faces_nearest(const std::string & fname, Regular_triangulation *rt)
{
	std::cout << "write " << fname << std::endl;

	std::ofstream out_file;
	out_file.open(fname);
	if (out_file.is_open())
	{
		int cnt = 0;
		for (auto it = rt->finite_faces_begin();
			it != rt->finite_faces_end(); ++it)
		{
			out_file
				<< it->estimated_cc() << " "
				<< it->nearest()
				<< std::endl;
			//std::cout << ++cnt << std::endl;
		}
		out_file.close();
		return true;
	}
	return false;
}

bool CReader::write_pts(const std::string & fname, std::set<Point_2> pts)
{
	std::cout << "write " << fname << std::endl;

	std::ofstream out_file;
	out_file.open(fname);
	if (out_file.is_open())
	{
		for (auto it = pts.begin(); it != pts.end(); ++it)
		{
			out_file << it->x() << " " << it->y() << std::endl;
		}
		out_file.close();
		return true;
	}
	return false;
}
