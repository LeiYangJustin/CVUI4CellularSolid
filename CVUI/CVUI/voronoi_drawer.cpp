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
	std::vector<std::vector<Point_2>> fitting_base_pts;
	pVD_->GetFittingBasePts(fitting_base_pts);

	cv_fitting_base_pts_list_.clear();
	for (int i = 0; i < fitting_base_pts.size(); i++)
	{
		std::vector<cv::Point> tmplist;
		for (int j = 1; j < fitting_base_pts[i].size(); j++)
		{
			cv::Point cvp1 = convert_to_cvPoint(fitting_base_pts[i][j]);
			tmplist.push_back(cvp1);
		}
		cv_fitting_base_pts_list_.push_back(tmplist);
	}
	

	//cv_fitting_base_pts_.clear();
	//for (int i = 0; i < fitting_base_pts.size(); i++)
	//{
	//	if (fitting_base_pts[i].size() == 1) {
	//		cv::Point cvp1 = convert_to_cvPoint(fitting_base_pts[i].front());
	//		cv::Point cvp2 = convert_to_cvPoint(fitting_base_pts[i].front());
	//		cv_fitting_base_pts_.push_back(std::make_pair(cvp1, cvp2));
	//	}
	//	else {
	//		for (int j = 1; j < fitting_base_pts[i].size(); j++)
	//		{
	//			cv::Point cvp1 = convert_to_cvPoint(fitting_base_pts[i][j]);
	//			cv::Point cvp2 = convert_to_cvPoint(fitting_base_pts[i][j-1]);
	//			cv_fitting_base_pts_.push_back(std::make_pair(cvp1, cvp2));
	//		}
	//	}
	//}
}

void CVoronoiDrawer::SetVD(CVoronoiDiagram* pVD)
{
	pVD_ = pVD;
	prepare_fitted_edges();
	prepare_voronoi_edges();
	prepare_fitting_base_pts();
}

void CVoronoiDrawer::SetImgData(CImgData* p_img_data)
{
	p_img_data_ = p_img_data;
}

cv::Point CVoronoiDrawer::convert_to_cvPoint(Point_2 p)
{
	return cv::Point(p.x(), p.y());
}
