// this is adapted from src code of de Goes et al. (2012)
// Blue Noise through Optimal Transport
// from http://fernandodegoes.org/

#include "scene.h"
#include "random.h"

void Scene::generate_random_sites_in_rectangle_box(const unsigned nb, const unsigned nx, const unsigned ny)
{
	std::vector<Point_2> points;
	double dx = m_domain.width();
	double dy = m_domain.height();
	while (points.size() != nb)
	{
		double x = random_double(-dx, dx);
		double y = random_double(-dy, dy);
		Point_2 p(x, y);
		if (m_domain.is_inside(p))
			points.push_back(p);
	}

	FT stepx = 2.0 * dx / nx;
	FT stepy = 2.0 * dy / ny;

	for (double bdy = -dy-stepy; bdy <= dy+stepy; bdy+=stepy)
	{
		Point_2 pmin(-dx - stepx, bdy);
		points.push_back(pmin);
		Point_2 pmax(dx + stepx, bdy);
		points.push_back(pmax);
	}
	for (double bdx = -dx - stepx; bdx <= dx + stepx; bdx += stepx)
	{
		Point_2 pmin(bdx, -dy - stepy);
		points.push_back(pmin);
		Point_2 pmax(bdx, dy + stepy);
		points.push_back(pmax);
	}

	std::vector<FT> weights(points.size(), 0.0);
	construct_triangulation(points, weights);
}

void Scene::generate_random_sites(const unsigned nb)
{
    std::vector<Point_2> points;
    double dx = m_domain.width();
    double dy = m_domain.height();
    while (points.size() != nb)
    {
        double x = random_double(-dx, dx);
        double y = random_double(-dy, dy);
        Point_2 p(x, y);
        if (m_domain.is_inside(p))
            points.push_back(p);
    }

    std::vector<FT> weights(points.size(), 0.0);
    construct_triangulation(points, weights);
}

void Scene::generate_regular_grid(const unsigned nx, const unsigned ny)
{
    double dx = m_domain.width();
    double dy = m_domain.height();
    FT stepx = 2.0 * dx / nx;
    FT stepy = 2.0 * dy / ny;
    std::vector<Point_2> points;
    for (unsigned i = 0; i < nx; ++i)
    {
        FT x = (i + 0.5)*stepx - dx;
        for (unsigned j = 0; j < ny; ++j)
        {
            FT y = (j + 0.5)*stepy - dy;
            Point_2 p(x, y);
            if (m_domain.is_inside(p))
                points.push_back(p);
        }
    }    
    std::vector<FT> weights(points.size(), 0.0);
    construct_triangulation(points, weights);
}

void Scene::generate_hextille_grid(const unsigned nb)
{
    double dx = m_domain.width();
    double dy = m_domain.height();
    FT step = 2.0 * dx / nb;
    std::vector<Point_2> points;
    for (unsigned i = 0; i < nb; ++i)
    {
        FT y = (i + 0.5)*step - dy;
        for (unsigned j = 0; j < nb; ++j)
        {
            FT x = (j + 0.25)*step - dx;
            if (i % 2 == 1) x += 0.5*step;
            Point_2 p(x, y);
            if (m_domain.is_inside(p))
                points.push_back(p);
        }
    }
    std::vector<FT> weights(points.size(), 0.0);
    construct_triangulation(points, weights);    
}
