#ifndef C_SKELETON_EXTRACTOR_H
#define C_SKELETON_EXTRACTOR_H

#include "prereq.h"
#include <vector>
#include "datatypedef.h"

// INCLUDE OPENCV DIR
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"

class CSkeletonExtractor
{
public:
	CSkeletonExtractor(cv::Mat gray_scale_src);
	~CSkeletonExtractor();

	std::vector<iPoint2> getSolidSkeleton();
	std::vector<iPoint2> getVoidSkeleton();

	void getVoidSkeletonSamples(std::vector<iPoint2> & X);

private:
	cv::Mat src_img_;
	cv::Mat rvs_img_;

	std::vector<iPoint2> solid_skeleton_;
	std::vector<iPoint2> void_skeleton_;

	// extraction of skeleton from a given image
	void extract_solid_morphological_skeleton();
	void extract_void_morphological_skeleton();
};

#endif