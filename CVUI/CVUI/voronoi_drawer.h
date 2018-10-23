#ifndef C_VORONOI_DRAWER_H
#define C_VORONOI_DRAWER_H

// INCLUDE OPENCV DIR
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include <opencv2/opencv.hpp>

#include <vector>
#include "../DataColle/voronoi_diagram.h"
#include "../DataColle/img_data.h"

class CVoronoiDrawer
{
public:
	CVoronoiDrawer();
	~CVoronoiDrawer();
	void SetVD(CVoronoiDiagram* pVD);
	void SetImgData(CImgData* p_img_data);
	CImgData* GetImgData() const { return p_img_data_; };

	//
	std::vector<std::pair<cv::Point, cv::Point>> GetVoronoiEdges()
	{
		return cv_voronoi_edges_;
	}
	std::vector<std::pair<cv::Point, cv::Point>> GetFittedEdges()
	{
		return cv_fitted_edges_;
	}
	//std::vector<std::pair<cv::Point, cv::Point>> GetFittingBasePoints()
	//{
	//	return cv_fitting_base_pts_;
	//}
	std::vector<std::vector<cv::Point>> GetFittingBasePointsList()
	{
		return cv_fitting_base_pts_list_;
	}

private:
	CVoronoiDiagram* pVD_;
	CImgData* p_img_data_;
	cv::Point convert_to_cvPoint(Point_2);

	std::vector<std::pair<cv::Point, cv::Point>> cv_voronoi_edges_;
	std::vector<std::pair<cv::Point, cv::Point>> cv_fitted_edges_;
	//std::vector<std::pair<cv::Point, cv::Point>> cv_fitting_base_pts_;
	std::vector<std::vector<cv::Point>> cv_fitting_base_pts_list_;
	//
	void prepare_voronoi_edges();
	void prepare_fitted_edges();
	void prepare_fitting_base_pts();
};

#endif