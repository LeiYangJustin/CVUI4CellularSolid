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

//void CVoronoiDiagram::UpdateTriangulation()
//{
//	typedef Regular_triangulation::Vertex_handle VertexHandle;
//	std::vector<WPoint> new_pts;
//	std::vector<WPoint> old_pts;
//	std::vector<bool> is_constrained_list;
//	for (auto viter = rt_->finite_vertices_begin(); viter != rt_->finite_vertices_end(); ++viter)
//	{
//		new_pts.push_back(viter->point());
//		//old_pts.push_back(WPoint(viter->info().point(), viter->info().weight()));
//		is_constrained_list.push_back(viter->info());
//	}
//
//	//
//	rt_->clear();
//	voronoi_segments_.clear();
//
//	for (int i = 0; i != new_pts.size(); ++i)
//	{
//		VertexHandle vh = rt_->insert(new_pts[i]);
//		vh->info() = is_constrained_list[i];
//	}
//	GetCroppedVoronoiSegments(voronoi_segments_);
//
//	//
//	std::cout << "num of hidden verts: " << rt_->number_of_hidden_vertices() << std::endl;
//}

void CVoronoiDiagram::SetSeedSites(const std::vector<WPoint>& wpts, const std::vector<bool> &is_constrained_list)
{
	rt_->clear();
	int cnt = 0;
	for (auto iter = wpts.begin(); iter != wpts.end(); ++iter, ++cnt)
	{
		auto vh = rt_->insert(*iter);
		vh->info() = is_constrained_list[cnt];
	}
}

void CVoronoiDiagram::GetCroppedVoronoiSegments(std::vector<Segment_2>& voronoi_segments)
{
	voronoi_segments.clear();
	Cropped_voronoi_from_delaunay vor(bbox_);
	//extract the cropped Voronoi diagram
	rt_->draw_dual(vor);
	voronoi_segments.insert(voronoi_segments.end(), vor.m_cropped_vd.begin(), vor.m_cropped_vd.end());
	//std::cout << vor.m_cropped_vd.size() << std::endl;
}

void CVoronoiDiagram::GetFittedSegments(std::vector<Segment_2>& fitting_segments)
{
	typedef Regular_triangulation::Finite_edges_iterator RT_Finite_Edge_Iter;
	typedef Regular_triangulation::Face_handle FaceHandle;

	std::vector<Segment_2> dual_edges;
	GetDualEdges(dual_edges);

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

void CVoronoiDiagram::GetCroppedTriangulatedSegments(std::vector<Segment_2>& triangulated_segments)
{
	triangulated_segments.clear();
	Cropped_voronoi_from_delaunay vor(bbox_);
	//extract the cropped Voronoi diagram
	rt_->draw_triangulation(vor);
	triangulated_segments.insert(triangulated_segments.end(), vor.m_cropped_vd.begin(), vor.m_cropped_vd.end());
}
