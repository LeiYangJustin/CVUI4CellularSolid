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
	clearMembers();
	rt_->clear();
	for (auto viter = wpts.begin(); viter != wpts.end(); ++viter)
	{
		rt_->insert(*viter);
	}
}

void CVoronoiDiagram::GetCroppedVoronoiSegments(Iso_rectangle_2 bbox, std::vector<Segment_2>& voronoi_segments)
{
	if (voronoi_segments_.size() == 0)
	{
		Cropped_voronoi_from_delaunay vor(bbox);
		//extract the cropped Voronoi diagram
		rt_->draw_dual(vor);
		voronoi_segments_.insert(voronoi_segments_.end(), vor.m_cropped_vd.begin(), vor.m_cropped_vd.end());
		std::cout << vor.m_cropped_vd.size() << std::endl;
	}
	voronoi_segments = voronoi_segments_;
}

void CVoronoiDiagram::GetFittedSegments(std::vector<Segment_2>& fitting_segments)
{
	typedef Regular_triangulation::Finite_edges_iterator RT_Finite_Edge_Iter;
	typedef Regular_triangulation::Face_handle FaceHandle;

	fitting_segments.clear();
	for (RT_Finite_Edge_Iter eit = rt_->finite_edges_begin(); eit != rt_->finite_edges_end(); ++eit)
	{
		int vind = eit->second;
		FaceHandle fh = eit->first;
		FaceHandle fh_oppo = fh->neighbor(vind);

		if (!rt_->is_infinite(fh) && !rt_->is_infinite(fh_oppo))
		{
			std::vector<Point_2> fh_plist = fh->info();
			Vector_2 fh_average_point(0, 0);
			for (int i = 0; i < fh_plist.size(); i++)
			{
				fh_average_point += (fh_plist[i] - Point_2(0, 0));
			}
			fh_average_point /= double(fh_plist.size());

			std::vector<Point_2> fh_oppo_plist = fh_oppo->info();
			Vector_2 fh_oppo_average_point(0, 0);
			for (int i = 0; i < fh_oppo_plist.size(); i++)
			{
				fh_oppo_average_point += (fh_oppo_plist[i] - Point_2(0, 0));
			}
			fh_oppo_average_point /= double(fh_oppo_plist.size());

			fitting_segments.push_back(Segment_2(Point_2(0,0)+fh_average_point, Point_2(0, 0)+fh_oppo_average_point));
		}
	}
}

void CVoronoiDiagram::GetFittingBasePts(std::vector<std::vector<Point_2>>& fitting_base_pts)
{
	typedef Regular_triangulation::Finite_faces_iterator RT_Finite_Face_Iter;
	typedef Regular_triangulation::Face_handle FaceHandle;
	
	fitting_base_pts.clear();
	for (RT_Finite_Face_Iter fit = rt_->finite_faces_begin(); fit != rt_->finite_faces_end(); ++fit)
	{
		if (!rt_->is_infinite(fit))
		{
			fitting_base_pts.push_back(fit->info());
		}
	}
}
