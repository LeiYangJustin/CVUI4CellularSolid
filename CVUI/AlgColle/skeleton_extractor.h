#ifndef C_SKELETON_EXTRACTOR_H
#define C_SKELETON_EXTRACTOR_H

#include "algprereq.h"
#include <vector>

// INCLUDE OPENCV DIR
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include <opencv2/opencv.hpp>

#include "../DataColle/img_data.h"

class ALGCOLLE_CLASS CSkeletonExtractor
{
public:
	CSkeletonExtractor(CImgData* img_data);
	~CSkeletonExtractor();

	std::vector<cv::Point> getSolidSkeleton();
	std::vector<cv::Point> getVoidSkeleton();

	void getVoidSkeletonSamples(std::vector<cv::Point> & X);
	void getImgData(CImgData* img_data) { img_data = img_data_; }
private:
	CImgData* img_data_;

	// extraction of skeleton from a given image
	void extract_morphological_skeleton(bool is_solid);
	void get_CVPoints_from_bwImg(cv::Mat bwImg, std::vector<cv::Point> &pts, bool use_boundary = true);

	// sampling
	int flood_fill_sampling(cv::Mat binaryGuide, std::vector<cv::Point> &samples, int radius = 10);
	//find joint and end points
	void find_critical_pts(cv::Mat binaryGuide, std::vector<cv::Point> skeleton,
		std::vector<cv::Point> &jointList, std::vector<cv::Point> &EndptList);
	void get_neighbors(cv::Point Apt, std::vector<cv::Point> &nns);
	/*
	a square element is like the following in the scanline ordering
	| 1 1 1 1 1 |
	| 1 1 1 1 1 |
	| 1 1 c 1 1 |
	| 1 1 1 1 1 |
	| 1 1 1 1 1 |
	*/
	void square_element(cv::Point Apt, int Rsize, std::vector<cv::Point> &nns);
	/*
	a outer skirt element is like the following without specific ordering
	| 1 1 1 1 1 |
	| 1 0 0 0 1 |
	| 1 0 c 0 1 |
	| 1 0 0 0 1 |
	| 1 1 1 1 1 |
	*/
	void outer_skirt_element(cv::Point Apt, int Rsize, std::vector<cv::Point> &nns);
	// skeleton thinning
	void skeleton_thinning(cv::Mat binaryInput, cv::Mat &binaryOutput);
	void skeleton_thinning_ZS_iter(cv::Mat& im, int iter);

	struct cvpt_compare {
		bool operator() (const cv::Point& lhs, const cv::Point& rhs) const {
			if (lhs.x != rhs.x)
				return lhs.x < rhs.x;
			else
				return lhs.y < rhs.y;
		}
	};
};

#endif