#ifndef C_SKELETON_EXTRACTOR_H
#define C_SKELETON_EXTRACTOR_H

#include "prereq.h"
#include <vector>

// INCLUDE OPENCV DIR
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"

class CSkeletonExtractor
{
public:
	CSkeletonExtractor(cv::Mat gray_scale_src);
	~CSkeletonExtractor();

	std::vector<cv::Point> getSolidSkeleton();
	std::vector<cv::Point> getVoidSkeleton();

	void getVoidSkeletonSamples(std::vector<cv::Point> & X);
	void getBoundingDomain(int & width, int & height);

private:
	cv::Mat src_img_;
	cv::Mat rvs_img_;

	std::vector<cv::Point> solid_skeleton_;
	std::vector<cv::Point> void_skeleton_;

	// extraction of skeleton from a given image
	void extract_solid_morphological_skeleton();
	void extract_void_morphological_skeleton();
};

#endif