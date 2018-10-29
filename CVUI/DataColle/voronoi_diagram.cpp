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


	//std::cout << "Testing GetCroppedVoronoiSegments: " << std::endl;
	//typedef Regular_triangulation::Face_handle FaceHandle;
	//typedef Regular_triangulation::Vertex_handle VertexHandle;
	//for (auto viter = voronoi_segments.begin(); viter != voronoi_segments.end(); ++viter)
	//{
	//	//// RT::locate is not a good way to find dual of the face 
	//	//// as the dual of a face may not lie in the face it corresponds to
	//	//FaceHandle src_fh = rt_->locate(WPoint(viter->source()));
	//	//FaceHandle tgt_fh = rt_->locate(WPoint(viter->target()));

	//	FaceHandle src_fh =NULL, tgt_fh=NULL;
	//	for (auto fiter = rt_->finite_faces_begin(); fiter != rt_->finite_faces_end(); ++fiter)
	//	{
	//		if (rt_->dual(fiter) == viter->source())
	//		{
	//			src_fh = fiter;
	//		}
	//		if (rt_->dual(fiter) == viter->target())
	//		{
	//			tgt_fh = fiter;
	//		}
	//		if (src_fh != NULL && tgt_fh != NULL)
	//		{
	//			if (!rt_->is_infinite(src_fh) && !rt_->is_infinite(tgt_fh)) {
	//				std::cout << "voronoi segment: " << viter->source() << " -> " << viter->target() << std::endl;
	//				std::cout << "and its corresponding face dual edge: "
	//					<< rt_->dual(src_fh) << " -> " << rt_->dual(tgt_fh) << std::endl;
	//			}
	//			break;
	//		}
	//	}
	//}

}

void CVoronoiDiagram::GetFittedSegments(std::vector<Segment_2>& fitting_segments)
{
	typedef Regular_triangulation::Finite_edges_iterator RT_Finite_Edge_Iter;
	typedef Regular_triangulation::Face_handle FaceHandle;

	std::vector<Segment_2> dual_edges;
	GetDualEdges(dual_edges);

	//fitting_segments.clear();
	//for (RT_Finite_Edge_Iter eit = rt_->finite_edges_begin(); eit != rt_->finite_edges_end(); ++eit)
	//{
	//	int vind = eit->second;
	//	FaceHandle fh = eit->first;
	//	FaceHandle fh_oppo = fh->neighbor(vind);

	//	if (!rt_->is_infinite(fh) && !rt_->is_infinite(fh_oppo))
	//	{
	//		std::vector<Point_2> fh_plist = fh->info();
	//		Point_2 fh_average_point(0, 0);
	//		for (int i = 0; i < fh_plist.size(); i++)
	//		{
	//			fh_average_point += (fh_plist[i] - Point_2(0, 0));
	//		}
	//		fh_average_point = Point_2(0, 0)+(fh_average_point-Point_2(0,0))/double(fh_plist.size());
	//		std::vector<Point_2> fh_oppo_plist = fh_oppo->info();
	//		Point_2 fh_oppo_average_point(0, 0);
	//		for (int i = 0; i < fh_oppo_plist.size(); i++)
	//		{
	//			fh_oppo_average_point += (fh_oppo_plist[i] - Point_2(0, 0));
	//		}
	//		fh_oppo_average_point = Point_2(0, 0) + (fh_oppo_average_point - Point_2(0, 0)) / double(fh_oppo_plist.size());
	//		fitting_segments.push_back(Segment_2(fh_average_point, fh_oppo_average_point));
	//	}
	//}

	fitting_segments.clear();
	for (auto deit = dual_edges.begin(); deit != dual_edges.end(); ++deit)
	{
		FaceHandle src_fh = NULL, tgt_fh = NULL;
		double src_dist = 1000, tgt_dist = 1000;
		for (auto fiter = rt_->finite_faces_begin(); fiter != rt_->finite_faces_end(); ++fiter)
		{
			Point_2 fPt = rt_->dual(fiter);
			if ((fPt-deit->source()).squared_length() < src_dist)
			{
				src_dist = (fPt - deit->source()).squared_length();
				src_fh = fiter;
			}
			if ((fPt- deit->target()).squared_length() < tgt_dist)
			{
				tgt_dist = (fPt - deit->target()).squared_length();
				tgt_fh = fiter;
			}
		}
		if (src_fh != NULL && tgt_fh != NULL)
		{
			std::vector<Point_2> src_plist = src_fh->info();
			Point_2 src_pt(0, 0);
			for (int i = 0; i < src_plist.size(); i++)
			{
				src_pt += (src_plist[i] - Point_2(0, 0));
			}
			src_pt = Point_2(0, 0)+(src_pt -Point_2(0,0))/double(src_plist.size());
			std::vector<Point_2> tgt_plist = tgt_fh->info();
			Point_2 tgt_pt(0, 0);
			for (int i = 0; i < tgt_plist.size(); i++)
			{
				tgt_pt += (tgt_plist[i] - Point_2(0, 0));
			}
			tgt_pt = Point_2(0, 0) + (tgt_pt - Point_2(0, 0)) / double(tgt_plist.size());
			fitting_segments.push_back(Segment_2(src_pt, tgt_pt));
		}
	}
}

void CVoronoiDiagram::GetDualEdges(std::vector<Segment_2>& dual_edges)
{
	typedef Regular_triangulation::Finite_edges_iterator RT_Finite_Edge_Iter;
	typedef Regular_triangulation::Face_handle FaceHandle;

	dual_edges.clear();
	for (RT_Finite_Edge_Iter eit = rt_->finite_edges_begin(); eit != rt_->finite_edges_end(); ++eit)
	{
		int vind = eit->second;
		FaceHandle fh = eit->first;
		FaceHandle fh_oppo = fh->neighbor(vind);

		if (!rt_->is_infinite(fh) && !rt_->is_infinite(fh_oppo))
		{
			dual_edges.push_back(Segment_2(rt_->dual(fh), rt_->dual(fh_oppo)));
		}
	}
}

void CVoronoiDiagram::GetFittingBasePtsMap(std::map<Point_2, std::vector<Point_2>>& fitting_base_pts_map)
{
	typedef Regular_triangulation::Finite_faces_iterator RT_Finite_Face_Iter;
	typedef Regular_triangulation::Face_handle FaceHandle;
	
	fitting_base_pts_map.clear();
	for (RT_Finite_Face_Iter fit = rt_->finite_faces_begin(); fit != rt_->finite_faces_end(); ++fit)
	{
		if (!rt_->is_infinite(fit))
		{
			fitting_base_pts_map.insert(std::make_pair(rt_->dual(fit), fit->info()));
		}
	}
}
