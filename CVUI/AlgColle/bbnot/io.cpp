// this is adapted from src code of de Goes et al. (2012)
// Blue Noise through Optimal Transport
// from http://fernandodegoes.org/

#include "scene.h"

#include <string>
#include <fstream>
#include <iostream>

void Scene::init_domain(unsigned nb)
{
    clear();
    if (nb != 4) 
    {
        std::cout << "only rectangular domains supported" << std::endl;
        return;
    }
    m_domain.init_rectangle(0.5, 0.5);
    m_domain.init_area();
}

//unsigned Scene::load(const QString& filename)
//{
//    std::vector<FT> weights;
//    std::vector<Point> points;
//
//    if (filename.contains(".dat", Qt::CaseInsensitive))
//    {
//        load_dat_file(filename, points);
//        weights.resize(points.size(), 0.0);
//    }
//    
//    if (filename.contains(".txt", Qt::CaseInsensitive))
//    {
//        load_txt_file(filename, points);
//        weights.resize(points.size(), 0.0);
//    }
//
//    if (!points.empty())
//    {
//        clear();
//        construct_triangulation(points, weights);
//        return m_vertices.size();
//    }
//    
//    std::cout << red << "try (.dat, .txt) file format" << white << std::endl;
//    return 0;
//}
//
//void Scene::load_dat_file(const QString& filename, std::vector<Point>& points) const
//{
//    std::ifstream ifs(qPrintable(filename));
//    Point point;
//    while (ifs >> point) points.push_back(point);
//    ifs.close();
//}
//
//void Scene::load_txt_file(const QString& filename, std::vector<Point>& points) const
//{
//    std::ifstream ifs(qPrintable(filename));
//
//    int numPoints;
//    ifs >> numPoints;
//
//    Point point;
//    while (ifs >> point) points.push_back(point);   
//    ifs.close();
//}
//
//void Scene::save(const QString& filename) const
//{
//    std::vector<Point> points;
//    collect_visible_points(points);
//    
//    if (filename.contains(".dat", Qt::CaseInsensitive))
//    {
//        save_dat(filename, points);
//        return;
//    }
//    
//    if (filename.contains(".txt", Qt::CaseInsensitive))
//    {
//        save_txt(filename, points);
//        return;
//    }
//    
//    if (filename.contains(".eps", Qt::CaseInsensitive))
//    {
//        save_eps(filename);
//        return;
//    }
//    
//    std::cout << red << "try (.dat, .txt, .eps) file format" << white << std::endl;
//}
//
//void Scene::save_dat(const QString& filename, const std::vector<Point>& points) const
//{
//    std::ofstream ofs(qPrintable(filename));
//    for (unsigned i = 0; i < points.size(); ++i)
//    {
//        ofs << points[i] << std::endl;
//    }
//    ofs.close();
//}
//
//void Scene::save_txt(const QString& filename, const std::vector<Point>& points) const
//{
//    std::ofstream ofs(qPrintable(filename));
//    ofs << points.size() << std::endl;
//    for (unsigned i = 0; i < points.size(); ++i)
//    {
//        ofs << points[i] << std::endl;
//    }
//    ofs.close();
//}
//
//void Scene::save_eps(const QString& filename) const
//{
//    double dx = m_domain.width();
//    double dy = m_domain.height();
//
//    double scale = 512.0;
//    double radius = 80.0 / m_vertices.size();
//
//    double min_x = -dx * scale;
//    double max_x =  dx * scale;
//    double min_y = -dy * scale;
//    double max_y =  dy * scale;
//
//
//    std::ofstream ofs(qPrintable(filename));
//    ofs << "%!PS-Adobe-3.1 EPSF-3.0\n";
//    ofs << "%%HiResBoundingBox: " 
//    << min_x << " " << min_y << " " << max_x << " " << max_y << std::endl;
//    ofs << "%%BoundingBox: " 
//    << min_x << " " << min_y << " " << max_x << " " << max_y << std::endl;
//    ofs << "%%CropBox: " 
//    << min_x << " " << min_y << " " << max_x << " " << max_y << "\n";
//    
//    ofs << "/radius { " << radius << " } def\n";
//    ofs << "/p { radius 0 360 arc closepath fill stroke } def\n";
//    ofs << "gsave " << scale << " " << scale << " scale\n";
//    ofs << "0 0 0 setrgbcolor" << std::endl;
//
//    for (unsigned i = 0; i < m_vertices.size(); ++i)
//    {
//        Vertex_handle vi = m_vertices[i];
//        if (vi->is_hidden()) continue;
//
//        const Point& pi = vi->get_position();
//        ofs << pi.x() << " " << pi.y() << " p" << std::endl;
//    }
//    ofs << "grestore" << std::endl;
//    ofs.close();
//}

// ADDED FUNCTIONS FOR VISUALIZATION
//
bool Scene::write_updated_triangulation(std::string fname, RT &rt)
{
	std::cout << "write updated triangulation" << std::endl;
	std::ofstream out_mesh_file;
	std::string foldername = "results\\";
	out_mesh_file.open(foldername + fname);
	
	if (out_mesh_file.is_open())
	{
		//std::map<Regular_triangulation::Vertex_handle, int> vh_id_map;
		//int cnt = 0;
		//for (auto vit = rt_->finite_vertices_begin(); vit != rt_->finite_vertices_end(); ++vit, ++cnt)
		//{
		//	vh_id_map.insert(std::make_pair(vit, cnt));
		//	out_mesh_file << cnt << " " << vit->point().point().x() << " " << vit->point().point().y() << std::endl;
		//}
		//cnt = 0;
		//for (auto fit = rt_->finite_faces_begin(); fit != rt_->finite_faces_end(); ++fit, ++cnt)
		//{
		//	int id0 = vh_id_map[fit->vertex(0)];
		//	int id1 = vh_id_map[fit->vertex(1)];
		//	int id2 = vh_id_map[fit->vertex(2)];
		//	out_mesh_file << cnt << " " << id0 << " " << id1 << " " << id2 << std::endl;
		//}
		//out_mesh_file.close();

		for (auto fit = rt.finite_faces_begin(); fit != rt.finite_faces_end(); ++fit)
		{
			out_mesh_file
				<< fit->vertex(0)->point().point() << " "
				<< fit->vertex(1)->point().point() << " "
				<< fit->vertex(2)->point().point() << " "
				<< std::endl;
		}
		out_mesh_file.close();
	}
	else {
		std::cout << "wrong at writing updated triangulation" << std::endl;
		return false;
	}

	std::cout << "save triangulation " << fname << std::endl;

	return true;
};