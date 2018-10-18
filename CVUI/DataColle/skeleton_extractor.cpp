#include "skeleton_extractor.h"

CSkeletonExtractor::CSkeletonExtractor(cv::Mat gray_scale_src)
{
	src_img_ = gray_scale_src;
}

CSkeletonExtractor::~CSkeletonExtractor()
{
}

std::vector<iPoint2> CSkeletonExtractor::getSolidSkeleton()
{
	return solid_skeleton_;
}

std::vector<iPoint2> CSkeletonExtractor::getVoidSkeleton()
{
	return void_skeleton_;
}

void CSkeletonExtractor::getVoidSkeletonSamples(std::vector<iPoint2>& X)
{
	// get void skeleton samples as the local extrema
}

void CSkeletonExtractor::extract_solid_morphological_skeleton()
{
	cv::Mat tmpImg = src_img_;

	// do extraction
}
void CSkeletonExtractor::extract_void_morphological_skeleton()
{
	cv::Mat tmpImg = rvs_img_;
	
	// do extraction
};
