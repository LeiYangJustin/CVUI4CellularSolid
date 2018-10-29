#include "voronoi_drawer.h"

CVoronoiDrawer::CVoronoiDrawer()
{
}

CVoronoiDrawer::~CVoronoiDrawer()
{
}

void CVoronoiDrawer::prepare_voronoi_edges()
{
	//
	int w, h;
	p_img_data_->GetBoundingDomain(w, h);

	//
	std::vector<Segment_2> voronoi_segments;
	pVD_->GetCroppedVoronoiSegments(Iso_rectangle_2(0, 0, w, h), voronoi_segments);

	//
	cv_voronoi_edges_.clear();
	for (int i = 0; i < voronoi_segments.size(); i++)
	{
		cv::Point cvp1 = convert_to_cvPoint(voronoi_segments[i].source());
		cv::Point cvp2 = convert_to_cvPoint(voronoi_segments[i].target());
		cv_voronoi_edges_.push_back(std::make_pair(cvp1, cvp2));
	}
}

void CVoronoiDrawer::prepare_fitted_edges()
{
	std::vector<Segment_2> fitted_segments;
	pVD_->GetFittedSegments(fitted_segments);

	cv_fitted_edges_.clear();
	for (int i = 0; i < fitted_segments.size(); i++)
	{
		cv::Point cvp1 = convert_to_cvPoint(fitted_segments[i].source());
		cv::Point cvp2 = convert_to_cvPoint(fitted_segments[i].target());
		cv_fitted_edges_.push_back(std::make_pair(cvp1, cvp2));
	}
}

void CVoronoiDrawer::prepare_fitting_base_pts()
{
	std::map<Point_2, std::vector<Point_2>> fitting_base_pts_map;
	pVD_->GetFittingBasePtsMap(fitting_base_pts_map);

	cv_fitting_base_pts_map_.clear();
	for (auto miter = fitting_base_pts_map.begin(); miter != fitting_base_pts_map.end(); ++miter)
	{
		std::vector<cv::Point> tmpVec;
		for (int i = 0; i < miter->second.size(); i++)
			tmpVec.push_back(convert_to_cvPoint(miter->second[i]));
		cv_fitting_base_pts_map_.insert(std::make_pair(convert_to_cvPoint(miter->first), tmpVec));
	}
}

void CVoronoiDrawer::prepare_dual_edges()
{
	std::vector<Segment_2> dual_edges;
	pVD_->GetDualEdges(dual_edges);

	cv_dual_edges_.clear();
	for (int i = 0; i < dual_edges.size(); i++)
	{
		cv::Point cvp1 = convert_to_cvPoint(dual_edges[i].source());
		cv::Point cvp2 = convert_to_cvPoint(dual_edges[i].target());
		cv_dual_edges_.push_back(std::make_pair(cvp1, cvp2));
	}
}

void CVoronoiDrawer::SetVD(CVoronoiDiagram* pVD)
{
	pVD_ = pVD;
	prepare_fitted_edges();
	prepare_voronoi_edges();
	prepare_dual_edges();
	//prepare_fitting_base_pts();
}

void CVoronoiDrawer::SetImgData(CImgData* p_img_data)
{
	p_img_data_ = p_img_data;
}

cv::Point CVoronoiDrawer::convert_to_cvPoint(Point_2 p)
{
	return cv::Point(p.x(), p.y());
}
