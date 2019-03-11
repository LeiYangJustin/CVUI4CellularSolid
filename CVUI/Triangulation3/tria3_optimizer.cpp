#include "tria3_optimizer.h"
#include "random3.h"

void CTria3Optimizer::generate_random_sites(const unsigned nb)
{
	std::vector<Point_3> points;
	double dw = domain_.width();
	while (points.size() != nb)
	{
		double x = random_double(-dw, dw);
		double y = random_double(-dw, dw);
		double z = random_double(-dw, dw);
		Point_3 p(x, y, z);
		if (domain_.is_inside(p))
			points.push_back(p);
	}

	std::vector<FT> weights(points.size(), 0.0);
	construct_triangulation(points, weights);
}

void CTria3Optimizer::lloyd_update()
{
	std::vector<Point_3> points;
	bool update = true;

	for (unsigned i = 0; i < vertices_.size(); ++i)
	{
		// update it to the barycenter of the corresponding voronoi cell
		Point_3 b = rt3_.compute_centroid(vertices_[i]);
		points.push_back(b);
	}

	update_positions(points);
	if (update)	update_triangulation();
}

void CTria3Optimizer::collect_sites(
	std::vector<Point_3> &points, 
	std::vector<FT> &weights) const
{
	for (unsigned i = 0; i < vertices_.size(); ++i)
	{
		Vertex_handle vi = vertices_[i];
		Point_3 pi = vi->get_position();
		points.push_back(pi);

		FT wi = 0.0;
		wi = vi->get_weight();
		weights.push_back(wi);
	}
}

void CTria3Optimizer::populate_vertices(std::vector<Point_3>& points, std::vector<FT>& weights)
{
	unsigned nb = 0;
	unsigned nsites = points.size();
	for (unsigned i = 0; i < nsites; ++i)
	{
		//std::cout << points[i] << std::endl;
		Vertex_handle vertex = insert_vertex(points[i], weights[i], nb);
		if (vertex == Vertex_handle()) continue;
		vertices_.push_back(vertex);
		nb++;
	}
}

void CTria3Optimizer::update_positions(std::vector<Point_3> points, bool clamp, bool hidden)
{
	unsigned j = 0;
	for (unsigned i = 0; i < vertices_.size(); ++i)
	{
		Vertex_handle vi = vertices_[i];
	
		//// skip hidden points
		//if (hidden && vi->is_hidden()) continue;

		// clamp
		Point_3 pi = points[j++];
		if (clamp) pi = domain_.clamp(pi);
		vi->set_position(pi);
	}
}

Vertex_handle CTria3Optimizer::insert_vertex(const Point_3 & point, const FT weight, const unsigned index)
{
	Weighted_point_3 wp(point, weight);
	Vertex_handle vertex = rt3_.insert(wp);

	if (vertex->get_index() != -1)
		return Vertex_handle();

	vertex->set_index(index);
	return vertex;
}

void CTria3Optimizer::construct_triangulation(std::vector<Point_3> points, std::vector<FT> weights)
{
	clear_triangulation();
	populate_vertices(points, weights);
	rt3_.pre_compute_cells();
	//compute_capacities(m_capacities);
}

void CTria3Optimizer::update_triangulation()
{
	std::vector<FT> weights;
	std::vector<Point_3> points;
	collect_sites(points, weights);
	construct_triangulation(points, weights);
}

void CTria3Optimizer::clear_triangulation()
{
	//m_capacities.clear();
	vertices_.clear();
	rt3_.clear();
}
