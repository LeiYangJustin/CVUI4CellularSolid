#ifndef C_IMG_DATA_H
#define C_IMG_DATA_H

#include "dataprereq.h"

#include <vector>
// INCLUDE OPENCV DIR
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include <opencv2/opencv.hpp>

class DATACOLLE_CLASS CImgData
{
public:
	CImgData(cv::Mat src_img);
	CImgData();
	~CImgData();
	bool ReadImgFromFile(const char* filename);
	//void SetSrcImg(cv::Mat src_img);

	cv::Mat GetVoidImg() { return void_img_; };
	cv::Mat GetSolidImg() { return solid_img_; };
	cv::Mat GetSrcImg() { return src_img_; };
	cv::Mat GetBackgroundImage() { make_background(); return background_img_; };

	void GetSolidSkeletonImg(cv::Mat &img) { convert_skeleton_to_img(true, img); };
	void GetVoidSkeletonImg(cv::Mat &img) { convert_skeleton_to_img(false, img); };

	std::vector<cv::Point>& GetSolidSkeleton() { return solid_skeletal_pts_; };
	std::vector<cv::Point>& GetVoidSkeleton() { return void_skeletal_pts_; };
	void SetSolidSkeleton(std::vector<cv::Point> solid_skeleton) { solid_skeletal_pts_ = solid_skeleton; };
	void SetVoidSkeletonImg(std::vector<cv::Point> void_skeleton) 
	{ has_void_skeleton_ = true;  void_skeletal_pts_ = void_skeleton; };

	bool HasVoidSkeleton() { return has_void_skeleton_; };
	void GetBoundingDomain(int & width, int & height);


	void find_solid_contours(std::vector<std::vector<double>> &edge_pts) { find_contours(edge_pts, solid_img_); };
	void find_void_contours(std::vector<std::vector<double>> &edge_pts) { find_contours(edge_pts, void_img_); };
	void get_two_distance_transform_fields(int &rows, int &cols,
		std::vector<double> &Sfield, std::vector<double> &Vfield);

private:
	cv::Mat src_img_;
	cv::Mat solid_img_;
	cv::Mat void_img_;	
	cv::Mat background_img_;
	bool has_void_skeleton_;

	std::vector<cv::Point> void_skeletal_pts_;
	std::vector<cv::Point> solid_skeletal_pts_;

	void binarize_src_img();
	void make_background();
	void convert_skeleton_to_img(bool is_solid, cv::Mat &img);
	void find_contours(std::vector<std::vector<double>> &edge_pts, const cv::Mat bw);
};

#endif // !C_IMG_DATA_H




