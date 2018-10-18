#include "voronoi_diagram.h"
#include <iostream>

CVoronoiDiagram::CVoronoiDiagram()
{
	rt_ = new Regular_triangulation;
}

CVoronoiDiagram::~CVoronoiDiagram()
{
	delete rt_;
}

void CVoronoiDiagram::SetSites(const std::vector<iPoint2>& pts)
{
	std::vector<WPoint> wps;
	for (int i = 0; i < pts.size(); i++)
	{
		BPoint bp(pts[i].x, pts[i].y);
		wps.push_back(WPoint(bp, 0.0));
	}
	rt_->insert(wps.begin(), wps.end());
}

void CVoronoiDiagram::addSite(const iPoint2 &p)
{
	BPoint bp(p.x, p.y);
	WPoint wp(bp, 0.0);
	rt_->insert(wp);
}

WPoint CVoronoiDiagram::pickSite(const iPoint2 &q) const
{
	// this can be replaced by the ANN search
	// find the nearest and then compute the distance
	typedef Regular_triangulation::Finite_vertices_iterator FVIT;
	FVIT fvit = rt_->finite_vertices_begin();
	for (; fvit != rt_->finite_vertices_end(); ++fvit) {
		if (abs(fvit->point().x() - q.x) < 5 && abs(fvit->point().y() - q.y) < 5)
		{
			// found
			break;
		}
	}
	return fvit->point();
}

void CVoronoiDiagram::setWeightToPickedSite(const BPoint &q, const double & w)
{
	// this can be replaced by the ANN search
	// find the nearest and then compute the distance
	typedef Regular_triangulation::Finite_vertices_iterator FVIT;
	FVIT fvit = rt_->finite_vertices_begin();
	for (; fvit != rt_->finite_vertices_end(); ++fvit) {
		if ( abs(fvit->point().x()-q.x()) < 5 && abs(fvit->point().y()-q.y()) < 5)
		{
			WPoint wp(fvit->point().point(), w);
			fvit->set_point(wp);
			break;
		}
	}
}

void CVoronoiDiagram::getVoronoiSegments(std::vector<std::pair<iPoint2, iPoint2>>& vseg_list)
{
	std::cout << "do something: get Voronoi segments" << std::endl;
}

void CVoronoiDiagram::getTriangulation(Regular_triangulation *rt)
{
	rt = rt_;
}
