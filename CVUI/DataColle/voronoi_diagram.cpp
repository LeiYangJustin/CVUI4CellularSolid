#include <iostream>

#include "voronoi_diagram.h"


CVoronoiDiagram::CVoronoiDiagram()
{
	rt_ = new Regular_triangulation;
}

CVoronoiDiagram::~CVoronoiDiagram()
{
	delete rt_;
}

Regular_triangulation* CVoronoiDiagram::getTriangulation()
{
	return rt_;
}

void CVoronoiDiagram::UpdateTriangulation(const std::vector<WPoint>& wpts)
{
	rt_->clear();
	for (auto viter = wpts.begin(); viter != wpts.end(); ++viter)
		rt_->insert(*viter);
}

void CVoronoiDiagram::GetCroppedVoronoiSegments(std::vector<std::pair<Point2, Point2>>& voronoi_segments)
{
	std::cout << "under construction" << std::endl;
}

//void CVoronoiDiagram::SetSites(const std::vector<BPoint>& pts)
//{
//	rt_->clear();
//	for (int i = 0; i < pts.size(); i++)
//	{
//		WPoint wp = WPoint(pts[i], 0.0);
//		rt_->insert(wp);
//	}
//}
//
//void CVoronoiDiagram::addSite(const BPoint &p)
//{
//	WPoint wp(p.x(), p.y());
//	rt_->insert(wp);
//}
//
//WPoint CVoronoiDiagram::pickSite(const BPoint &q) const
//{
//	// this can be replaced by the ANN search
//	// find the nearest and then compute the distance
//	typedef Regular_triangulation::Finite_vertices_iterator FVIT;
//	FVIT fvit = rt_->finite_vertices_begin();
//	for (; fvit != rt_->finite_vertices_end(); ++fvit) {
//		if (abs(fvit->point().x() - q.x()) < 5 && abs(fvit->point().y() - q.y()) < 5)
//		{
//			// found
//			break;
//		}
//	}
//	return fvit->point();
//}
//
//void CVoronoiDiagram::setWeightToPickedSite(const BPoint &q, const double & w)
//{
//	// this can be replaced by the ANN search
//	// find the nearest and then compute the distance
//	typedef Regular_triangulation::Finite_vertices_iterator FVIT;
//	FVIT fvit = rt_->finite_vertices_begin();
//	for (; fvit != rt_->finite_vertices_end(); ++fvit) {
//		if ( abs(fvit->point().x()-q.x()) < 5 && abs(fvit->point().y()-q.y()) < 5)
//		{
//			WPoint wp(fvit->point().point(), w);
//			fvit->set_point(wp);
//			break;
//		}
//	}
//}
//
//void CVoronoiDiagram::getVoronoiSegments(std::vector<std::pair<WPoint, WPoint>>& vseg_list)
//{
//	std::cout << "do something: get Voronoi segments" << std::endl;
//}

