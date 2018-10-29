#include "reconstructor.h"
#include <map>
#include <queue>
#include "../AlgColle/geo_calculator.h"
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
	//typedef Regular_triangulation::Finite_edges_iterator RT_Finite_Edge_Iter;
	typedef Regular_triangulation::Face_handle FaceHandle;
	typedef Regular_triangulation::Vertex_handle VertexHandle;

	Regular_triangulation* rt = pVD->getTriangulation();
	// get voronoi segments
	std::vector<Segment_2> voronoi_segments;
	pVD->GetCroppedVoronoiSegments(Iso_rectangle_2(0,0,width_,height_), voronoi_segments);
	//pVD->GetDualEdges(voronoi_segments);

	// Initialize the recon_edge_list with voronoi edges
	// assign label and distance to each reconstruction point
	// make point sets according to the labels
	// compute reconstruction_error_
	std::vector<CReconstructionEdge> recon_edge_list;
	for (int i = 0; i < voronoi_segments.size(); i++)
	{
		CReconstructionEdge reconEdge;
		reconEdge.set_id(i);
		reconEdge.set_voronoi_edge(voronoi_segments[i]);
		recon_edge_list.push_back(reconEdge);
	}
	std::map<int, std::vector<CReconstructionPoint>> map_lid_to_rps;
	for (auto rp_iter = reconPointList_.begin(); rp_iter != reconPointList_.end(); ++rp_iter)
	{
		// assign label and distance to each reconstruction point
		for (int lid = 0; lid < voronoi_segments.size(); ++lid)
		{
			double dist = rp_iter->DistToSegment(voronoi_segments[lid]);
			if (rp_iter->GetDist() > dist)
			{
				rp_iter->SetLabelDist(lid, dist);
			}
		}
		// add the point to the point set
		if (map_lid_to_rps.find(rp_iter->GetLabel()) != map_lid_to_rps.end())
		{
			map_lid_to_rps[rp_iter->GetLabel()].push_back(*rp_iter);
		}
		else {
			std::vector<CReconstructionPoint> emptyList;
			emptyList.push_back(*rp_iter);
			map_lid_to_rps.insert(std::make_pair(rp_iter->GetLabel(), emptyList));
		}
	}

	// Check if the initial assignment meets the good fit criteria
	std::queue<CReconstructionPoint> orphan_rpts_queue;
	for (auto miter = map_lid_to_rps.begin(); miter != map_lid_to_rps.end(); ++miter)
	{
		if (miter->first == -1)
		{
			std::vector<CReconstructionPoint> tmp_pts = miter->second;
			//std::cout << tmp_pts.size() << std::endl;
			for (auto iter = tmp_pts.begin(); iter != tmp_pts.end(); ++iter)
				orphan_rpts_queue.push(CReconstructionPoint(iter->point_, -1));
		}
		else 
		{
			CReconstructionEdge rc_edge = recon_edge_list[miter->first];
			rc_edge.set_point_list(miter->second);
			if (!rc_edge.is_good_fit() || miter->second.size() < 5)
			{
				std::vector<CReconstructionPoint> tmp_pts;
				rc_edge.get_point_list(tmp_pts);
				//std::cout << tmp_pts.size() << std::endl;
				for (auto iter = tmp_pts.begin(); iter != tmp_pts.end(); ++iter)
					orphan_rpts_queue.push(CReconstructionPoint(iter->point_, -1));
			}
			else
			{
				recon_edge_list[miter->first].set_point_list(miter->second);
			}
		}
	}

	// Iteratively to add orphan pts to the existing fitting model
	std::cout << "orphan_rpts_queue_size before: " << orphan_rpts_queue.size() << std::endl;
	int orphan_rpts_queue_size =  orphan_rpts_queue.size();
	int maxIter = 3 * orphan_rpts_queue_size;
	while(--maxIter && !orphan_rpts_queue.empty())
	{
		CReconstructionPoint rc_pt(orphan_rpts_queue.front().point_);
		orphan_rpts_queue.pop();
		for (auto viter = recon_edge_list.begin(); viter != recon_edge_list.end(); ++viter)
		{
			if (viter->can_add_point(rc_pt))
			{
				double dist = rc_pt.DistToSegment(viter->get_voronoi_edge());
				if (dist < rc_pt.GetDist())	{
					rc_pt.SetLabelDist(viter->get_id(), dist);
				}
			}
		}
		if (rc_pt.GetLabel() != -1)
		{
			int id = rc_pt.GetLabel();
			recon_edge_list[id].add_point(rc_pt);
		}
		else
			orphan_rpts_queue.push(rc_pt);
	}
	std::cout << "orphan_rpts_queue_size: " << orphan_rpts_queue.size() << "; iteration: " << maxIter << std::endl;

	// Update the face centroids
	for (auto viter = recon_edge_list.begin(); viter != recon_edge_list.end(); ++viter)
	{
		if (viter->is_good_fit())
		{
			Point_2 fitp1 = viter->get_fitted_edge().source();
			Point_2 fitp2 = viter->get_fitted_edge().target();
			Point_2 src_pt = viter->get_voronoi_edge().source();
			Point_2 tgt_pt = viter->get_voronoi_edge().target();
			//std::cerr << "still good" << std::endl;
			FaceHandle src_fh = NULL, tgt_fh = NULL;
			for (auto fiter = rt->finite_faces_begin(); fiter != rt->finite_faces_end(); ++fiter)
			{
				if (rt->dual(fiter) == src_pt) {
					src_fh = fiter;
				}
				if (rt->dual(fiter) == tgt_pt) {
					tgt_fh = fiter;
				}
				if (src_fh != NULL && tgt_fh != NULL) {
					if (CGAL::compare_distance_to_point(src_pt, fitp1, fitp2) == CGAL::SMALLER) {
						src_fh->info().push_back(fitp1);
						tgt_fh->info().push_back(fitp2);
					}
					else {
						src_fh->info().push_back(fitp2);
						tgt_fh->info().push_back(fitp1);
					}
					break;
				}
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

	// Update reconstruction error
	reconstruction_error_ = 0.0;
	for (auto viter = recon_edge_list.begin(); viter != recon_edge_list.end(); ++viter)
	{
		std::vector<CReconstructionPoint> tmp_pts;
		viter->get_point_list(tmp_pts);
		for (int i = 0; i < tmp_pts.size(); i++)
		{
			reconstruction_error_ += tmp_pts[i].DistToSegment(viter->get_fitted_edge());
		}
	}


	// Update reconstruction point for visualization
	map_lid_to_rps.clear();
	reconPointList_.clear();
	for (auto viter = recon_edge_list.begin(); viter != recon_edge_list.end(); ++viter)
	{
		std::vector<CReconstructionPoint> tmp_pts;
		viter->get_point_list(tmp_pts);
		if (tmp_pts.size() > 0) {
			reconPointList_.insert(reconPointList_.end(), tmp_pts.begin(), tmp_pts.end());
		}
	}
	std::cout << "reconPointList_ size: " << reconPointList_.size() << std::endl;
	std::queue<CReconstructionPoint> tmp_orphan_rpts_queue = orphan_rpts_queue;
	while (!tmp_orphan_rpts_queue.empty())
	{
		reconPointList_.push_back(CReconstructionPoint(tmp_orphan_rpts_queue.front().point_, -1));
		tmp_orphan_rpts_queue.pop();
	}
	std::cout << "reconPointList_ size: " << reconPointList_.size() << std::endl;
}

void CReconstructor::Update2(CVoronoiDiagram * pVD)
{
	//typedef Regular_triangulation::Finite_edges_iterator RT_Finite_Edge_Iter;
	typedef Regular_triangulation::Face_handle FaceHandle;
	typedef Regular_triangulation::Vertex_handle VertexHandle;

	// get voronoi segments
	Regular_triangulation* rt = pVD->getTriangulation();
	std::vector<Segment_2> voronoi_segments;
	pVD->GetCroppedVoronoiSegments(Iso_rectangle_2(0, 0, width_, height_), voronoi_segments);

	// STEP 2
	// assign label and distance to each reconstruction point
	// make point sets according to the labels
	// compute reconstruction_error_
	std::map<int, std::list<CReconstructionPoint>> map_lid_to_rplist;
	reconstruction_error_ = 0.0;
	for (auto rp_iter = reconPointList_.begin(); rp_iter != reconPointList_.end(); ++rp_iter)
	{
		// assign label and distance to each reconstruction point
		for (int lid = 0; lid < voronoi_segments.size(); ++lid)
		{
			double dist = rp_iter->DistToSegment(voronoi_segments[lid]);
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
	
	//fit_line_segment_to_labeled_pts
	// fit a line segment to each of the point sets
	for (auto miter = map_lid_to_rplist.begin(); miter != map_lid_to_rplist.end(); ++miter)
	{
		if (miter->second.size() > 20) {
			Point_2 fitp1, fitp2;
			if (update_edge_by_fitting_new(fitp1, fitp2, miter->second) < 0.3)
			{
				//std::cerr << "this is not a goodfit because of low linearity" << std::endl;
				continue;
			}
			Point_2 src_pt = voronoi_segments[miter->first].source();
			Point_2 tgt_pt = voronoi_segments[miter->first].target();
			if (fabs(CGeoCalculator::compute_cos_angle(Segment_2(fitp1, fitp2), voronoi_segments[miter->first])) < 0.7) // cosd(45)
			{
				//std::cerr << "this is not a goodfit because of dis-alignment" << std::endl;
				continue;
			}
			//std::cerr << "still good" << std::endl;
			FaceHandle src_fh = NULL, tgt_fh = NULL;
			for (auto fiter = rt->finite_faces_begin(); fiter != rt->finite_faces_end(); ++fiter)
			{
				if (rt->dual(fiter) == src_pt) {
					src_fh = fiter;
				}
				if (rt->dual(fiter) == tgt_pt) {
					tgt_fh = fiter;
				}
				if (src_fh != NULL && tgt_fh != NULL) {
					if (CGAL::compare_distance_to_point(src_pt, fitp1, fitp2) == CGAL::SMALLER) {
						src_fh->info().push_back(fitp1);
						tgt_fh->info().push_back(fitp2);
					}
					else {
						src_fh->info().push_back(fitp2);
						tgt_fh->info().push_back(fitp1);
					}
					break;
				}
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

double CReconstructor::GetReconstructionError()
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

double CReconstructor::update_edge_by_fitting_new(Point_2 & fitp1, Point_2 & fitp2, std::list<CReconstructionPoint> rp_list)
{
	std::vector<Point_2> points;
	for (auto liter = rp_list.begin(); liter != rp_list.end(); ++liter)
	{
		points.push_back(liter->point_);
	}

	// fit line 
	Line_2 line;
	double goodfit = linear_least_squares_fitting_2(points.begin(), points.end(), line, CGAL::Dimension_tag<0>());

	// finding the two ends
	// Projection
	Point_2 center = centroid(points.begin(), points.end());
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

	return goodfit;
}

void CReconstructor::fit_line_segment_to_labeled_pts(
	CVoronoiDiagram *pVD,
	std::vector<Segment_2> voronoi_segments,
	std::map<int, std::list<CReconstructionPoint>> map_lid_to_rplist)
{
	typedef Regular_triangulation::Face_handle FaceHandle;

	Regular_triangulation *rt = pVD->getTriangulation();

	// fit a line segment to each of the point sets
	for (auto miter = map_lid_to_rplist.begin(); miter != map_lid_to_rplist.end(); ++miter)
	{
		if (miter->second.size() > 20) {
			Point_2 fitp1, fitp2;
			if (update_edge_by_fitting_new(fitp1, fitp2, miter->second) < 0.3)
			{
				//std::cerr << "this is not a goodfit because of low linearity" << std::endl;
				continue;
			}
			Point_2 src_pt = voronoi_segments[miter->first].source();
			Point_2 tgt_pt = voronoi_segments[miter->first].target();
			if (fabs(CGeoCalculator::compute_cos_angle(Segment_2(fitp1, fitp2), voronoi_segments[miter->first])) < 0.7) // cosd(45)
			{
				//std::cerr << "this is not a goodfit because of dis-alignment" << std::endl;
				continue;
			}
			//std::cerr << "still good" << std::endl;
			FaceHandle src_fh = NULL, tgt_fh = NULL;
			for (auto fiter = rt->finite_faces_begin(); fiter != rt->finite_faces_end(); ++fiter)
			{
				if (rt->dual(fiter) == src_pt) {
					src_fh = fiter;
				}
				if (rt->dual(fiter) == tgt_pt) {
					tgt_fh = fiter;
				}
				if (src_fh != NULL && tgt_fh != NULL) {
					if (CGAL::compare_distance_to_point(src_pt, fitp1, fitp2) == CGAL::SMALLER) {
						src_fh->info().push_back(fitp1);
						tgt_fh->info().push_back(fitp2);
					}
					else {
						src_fh->info().push_back(fitp2);
						tgt_fh->info().push_back(fitp1);
					}
					break;
				}
			}
		}
	}

	//// Test
	//for (auto miter = map_lid_to_rplist.begin(); miter != map_lid_to_rplist.end(); ++miter)
	//{
	//	Point_2 fitp1, fitp2;
	//	if (update_edge_by_fitting_new(fitp1, fitp2, miter->second) < 0.3)
	//	{
	//		//std::cerr << "this is not a goodfit because of low linearity" << std::endl;
	//		continue;
	//	}
	//	Point_2 src_pt = voronoi_segments[miter->first].source();
	//	Point_2 tgt_pt = voronoi_segments[miter->first].target();

	//	Point_2 test_src_pt(388.076, 272.137);
	//	Point_2 test_tgt_pt(370.524, 240.543);

	//	if ((src_pt - test_src_pt).squared_length() < 0.05 && (tgt_pt - test_tgt_pt).squared_length())
	//	{
	//		std::vector<Point_2> points;
	//		for (auto liter = miter->second.begin(); liter != miter->second.end(); ++liter)
	//		{
	//			points.push_back(liter->point_);
	//		}
	//		// fit line 
	//		Line_2 line;
	//		double goodfit = linear_least_squares_fitting_2(points.begin(), points.end(), line, CGAL::Dimension_tag<0>());
	//		Point_2 center = centroid(points.begin(), points.end());
	//		Vector_2 vdir = line.to_vector();
	//		vdir /= sqrt(vdir.squared_length());

	//		double min_proj_para = 100, max_proj_para = -100;
	//		for (int i = 0; i < points.size(); i++)
	//		{
	//			Vector_2 v1 = points[i] - center;
	//			double tmp_proj_para;
	//			double sign = v1*vdir;
	//			tmp_proj_para = v1*vdir / sqrt(vdir.squared_length());
	//			if (min_proj_para > tmp_proj_para)
	//				min_proj_para = tmp_proj_para;
	//			if (max_proj_para < tmp_proj_para)
	//				max_proj_para = tmp_proj_para;
	//		}
	//		fitp1 = center + min_proj_para*vdir;
	//		fitp2 = center + max_proj_para*vdir;

	//		std::cout << "Voronoi edge: " << src_pt << " -> " << tgt_pt << std::endl;
	//		std::cout << " fitted edge: " << fitp1 << " -> " << fitp2 << std::endl;
	//		std::cout << "   paramters: " << min_proj_para << " -> " << max_proj_para << std::endl;
	//		std::cout << " center: " << center << " and vdir " << vdir << std::endl;
	//		std::cout << "assigned rc pts: ";
	//		std::copy(points.begin(), points.end(), std::ostream_iterator<Point_2>(std::cout, "/"));
	//		std::cout << std::endl;
	//		std::cout << std::endl;
	//	}
	//	else if ((src_pt - test_tgt_pt).squared_length() < 0.05 && (tgt_pt -  test_src_pt).squared_length())
	//	{
	//		std::vector<Point_2> points;
	//		for (auto liter = miter->second.begin(); liter != miter->second.end(); ++liter)
	//		{
	//			points.push_back(liter->point_);
	//		}		
	//		// fit line 
	//		Line_2 line;
	//		double goodfit = linear_least_squares_fitting_2(points.begin(), points.end(), line, CGAL::Dimension_tag<0>());
	//		Point_2 center = centroid(points.begin(), points.end());
	//		Vector_2 vdir = line.to_vector();
	//		vdir /= sqrt(vdir.squared_length());

	//		std::cout << "Voronoi edge: " << src_pt << " -> " << tgt_pt << std::endl;
	//		std::cout << " fitted edge: " << fitp1 << " -> " << fitp2 << std::endl;
	//		std::cout << " center: " << center << " and vdir " << vdir << std::endl;
	//		std::cout << "assigned rc pts: ";
	//		std::copy(points.begin(), points.end(), std::ostream_iterator<Point_2>(std::cout, "/"));
	//		std::cout << std::endl;
	//		std::cout << std::endl;
	//	}
	//}
}
