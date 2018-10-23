#include "reconstructor.h"
#include <map>
#include "geo_calculator.h"
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/centroid.h>
CReconstructor::CReconstructor() :
	reconstruction_error_(0.0)
{
}


CReconstructor::~CReconstructor()
{
}

void CReconstructor::SetBBox(int w, int h)
{
	width_ = w;
	height_ = h;
}

void CReconstructor::SetReconstructionPts(std::vector<Point_2> P)
{
	for (int i = 0; i < P.size(); i++)
	{
		CReconstructionPoint rp(P[i]);
		reconPointList_.push_back(rp);
	}
}

void CReconstructor::GetReconstructionPts(std::vector<Point_2>& RPts, std::vector<int>& RPt_labels)
{
	RPts.clear();
	RPt_labels.clear();
	for (int i = 0; i < reconPointList_.size(); i++)
	{
		RPts.push_back(reconPointList_[i].point_);
		RPt_labels.push_back(reconPointList_[i].GetLabel());
	}
}

void CReconstructor::Update(CVoronoiDiagram * pVD)
{
	//std::cout << "reconstruction" << std::endl;

	//typedef Regular_triangulation::Finite_edges_iterator RT_Finite_Edge_Iter;
	typedef Regular_triangulation::Face_handle FaceHandle;
	typedef Regular_triangulation::Vertex_handle VertexHandle;

	// get voronoi segments
	Regular_triangulation* rt = pVD->getTriangulation();
	std::vector<Segment_2> voronoi_segments;
	pVD->GetCroppedVoronoiSegments(Iso_rectangle_2(0, 0, width_, height_), voronoi_segments);

	// assign label and distance to each reconstruction point
	// make point sets according to the labels
	// compute reconstruction_error_
	reconstruction_error_ = 0.0;
	std::map<int, std::list<CReconstructionPoint>> map_lid_to_rplist;
	for (auto rp_iter = reconPointList_.begin(); rp_iter != reconPointList_.end(); ++rp_iter)
	{
		// assign label and distance to each reconstruction point
		for (int lid = 0; lid < voronoi_segments.size(); ++lid)
		{
			double dist = rp_iter->SqrDistToSegment(voronoi_segments[lid]);
			if (rp_iter->GetDist() > dist)
			{
				rp_iter->SetLabelDist(lid, dist);
			}
		}

		// add the point to the point set
		if (map_lid_to_rplist.find(rp_iter->GetLabel()) != map_lid_to_rplist.end())
		{
			map_lid_to_rplist[rp_iter->GetLabel()].push_back(*rp_iter);
		}
		else {
			std::list<CReconstructionPoint> emptyList;
			emptyList.push_back(*rp_iter);
			map_lid_to_rplist.insert(std::make_pair(rp_iter->GetLabel(), emptyList));
		}

		// compute reconstruction_error_
		reconstruction_error_ += rp_iter->GetDist();
	}

	// fit a line segment to each of the point sets
	for (auto miter = map_lid_to_rplist.begin(); miter != map_lid_to_rplist.end(); ++miter)
	{
		//std::cout << miter->second.size() << std::endl;
		if (miter->second.size() > 50) {
			Point_2 fitp1, fitp2;
			update_edge_by_fitting_new(fitp1, fitp2, miter->second);

			Point_2 src_pt = voronoi_segments[miter->first].source();
			Point_2 tgt_pt = voronoi_segments[miter->first].target();
			FaceHandle src_fh = rt->locate(WPoint(src_pt));
			FaceHandle tgt_fh = rt->locate(WPoint(tgt_pt));

			if (CGAL::compare_distance_to_point(src_pt, fitp1, fitp2) == CGAL::SMALLER)
			{
				src_fh->info().push_back(fitp1);
				tgt_fh->info().push_back(fitp2);
				std::cout << "fitted edge: (" << fitp1 << "), (" << fitp2 << ")" << std::endl;
				std::cout << "voronoi edge: (" << src_pt << "), (" << tgt_pt << ")\n" << std::endl;
			}
			else {
				src_fh->info().push_back(fitp2);
				tgt_fh->info().push_back(fitp1);
				std::cout << "fitted edge: (" << fitp2 << "), (" << fitp1 << ")" << std::endl;
				std::cout << "voronoi edge: (" << src_pt << "), (" << tgt_pt << ")\n" << std::endl;
			}
		}
	}

	// process those faces with no fitting information
	for (auto fit = rt->finite_faces_begin(); fit != rt->finite_faces_end(); ++fit)
	{
		if (fit->info().size() == 0)
		{
			fit->info().push_back(rt->dual(fit));
		}
	}
}

double CReconstructor::GetAccuracy()
{
	return reconstruction_error_;
}

void CReconstructor::update_edge_by_fitting(Point_2 &fitp1, Point_2 &fitp2,
	std::list<CReconstructionPoint> rp_list)
{
	// prepare data
	std::vector<Point_2> plist;
	Vector_2 cvec(0, 0);
	for (auto siter = rp_list.begin(); siter != rp_list.end(); ++siter)
	{
		plist.push_back(siter->point_);
		cvec +=  (siter->point_ - Point_2(0,0));
	}
	cvec /= plist.size();
	
	// fitting a line in the first place;
	// PCA
	std::vector<std::vector<double>> eigVecs;
	std::vector<double> eigVals;
	CGeoCalculator::compute_PCA_with_SVD(plist, eigVecs, eigVals);
	Vector_2 vdir(eigVecs[0][0], eigVecs[0][1]);

	// finding the two ends
	// Projection
	Point_2 center(cvec.x(), cvec.y());
	double min_proj_para = 100, max_proj_para = -100;
	for (int i = 0; i < plist.size(); i++)
	{
		double tmp_proj_para = CGeoCalculator::projection_of_point_to_line(plist[i], center, center + vdir);
		if (min_proj_para > tmp_proj_para)
			min_proj_para = tmp_proj_para;
		if (max_proj_para < tmp_proj_para)
			max_proj_para = tmp_proj_para;
	}
	fitp1 = center + min_proj_para*vdir;
	fitp2 = center + max_proj_para*vdir;
}

void CReconstructor::update_edge_by_fitting_new(Point_2 & fitp1, Point_2 & fitp2, std::list<CReconstructionPoint> rp_list)
{
	std::vector<Point_2> points;
	for (auto liter = rp_list.begin(); liter != rp_list.end(); ++liter)
	{
		points.push_back(liter->point_);
	}

	// fit line 
	Line_2 line;
	linear_least_squares_fitting_2(points.begin(), points.end(), line, CGAL::Dimension_tag<0>());

	Point_2 center = centroid(points.begin(), points.end());

	// finding the two ends
	// Projection
	Vector_2 vdir = line.to_vector();
	double min_proj_para = 100, max_proj_para = -100;
	for (int i = 0; i < points.size(); i++)
	{
		double tmp_proj_para = CGeoCalculator::projection_of_point_to_line(points[i], center, center + vdir);
		if (min_proj_para > tmp_proj_para)
			min_proj_para = tmp_proj_para;
		if (max_proj_para < tmp_proj_para)
			max_proj_para = tmp_proj_para;
	}
	fitp1 = center + min_proj_para*vdir;
	fitp2 = center + max_proj_para*vdir;
}
