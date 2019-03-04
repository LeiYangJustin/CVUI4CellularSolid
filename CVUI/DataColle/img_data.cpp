#include "img_data.h"



CImgData::CImgData(cv::Mat src_img)
{
	src_img_ = src_img;
	// BINARIZE THE EXEMPLAR
	binarize_src_img();
}

CImgData::CImgData() : has_void_skeleton_(false)
{
}


CImgData::~CImgData()
{
	void_skeletal_pts_.clear();
	solid_skeletal_pts_.clear();
}

bool CImgData::ReadImgFromFile(const char* filename)
{
	cv::Mat src_img = cv::imread(filename);
	if (src_img.data == NULL) {
		std::cerr << "Cannot successfully read the image file" << std::endl;
		return false;
	}

#ifdef _DEBUG
	int dRows = 101;
	//int dCols = 101;
	std::cout << "We are in the debugging mode" << std::endl;
#else 
	int dRows = 401;
	//int dCols = 501;
#endif
	// resize the image
	int dCols = double(dRows) / double(src_img_.rows)*double(src_img_.cols);
	cv::resize(src_img, src_img_, cv::Size(dCols, dRows));

	// BINARIZE THE EXEMPLAR
	binarize_src_img();
	return true;
}


void CImgData::GetBoundingDomain(int & width, int & height)
{
	width = src_img_.cols;
	height = src_img_.rows;
}

void CImgData::get_two_distance_transform_fields(int &rows, int &cols,
	std::vector<double> &Sfield, std::vector<double> &Vfield)
{
	//
	cols = solid_img_.cols;
	rows = solid_img_.rows;

	//
	cv::Mat solid_dist(solid_img_.size(), CV_32F);
	cv::distanceTransform(solid_img_, solid_dist, cv::DIST_L2, cv::DIST_MASK_3);
	Sfield.clear();
	for (int ix = 0; ix < solid_dist.cols; ix++)
	{
		for (int iy = 0; iy < solid_dist.rows; iy++)
		{
			Sfield.push_back(-solid_dist.at<float>(iy, ix));
		}
	}

	//
	cv::Mat void_dist(void_img_.size(), CV_32F);
	cv::distanceTransform(void_img_, void_dist, cv::DIST_L2, cv::DIST_MASK_5);
	Vfield.clear();
	for (int ix = 0; ix < void_dist.cols; ix++)
	{
		for (int iy = 0; iy < void_dist.rows; iy++)
		{
			Vfield.push_back(-void_dist.at<float>(iy, ix));
		}
	}
}

void CImgData::binarize_src_img()
{
	// BINARIZE THE EXEMPLAR
	cv::Mat grayExplImg(src_img_.size(), CV_8UC1);
	cv::cvtColor(src_img_, grayExplImg, cv::COLOR_BGR2GRAY);
	cv::blur(grayExplImg, solid_img_, cv::Size(3, 3));
	cv::threshold(solid_img_, solid_img_, 50, 255, cv::THRESH_BINARY | cv::THRESH_OTSU);
	cv::bitwise_not(solid_img_, void_img_);
	//
	background_img_ = solid_img_;
}

void CImgData::make_background()
{
	std::vector<cv::Mat> channels(3);
	cv::Mat chnl2 = this->GetSolidImg().clone();
	cv::Mat chnl1 = this->GetSolidImg().clone();
	cv::Mat chnl0 = this->GetSolidImg().clone();

	////
	//cv::subtract(chnl0, skelImg_secondary_, chnl0);
	//cv::subtract(chnl2, skelImg_secondary_, chnl2);
	//cv::add(chnl1, skelImg_secondary_, chnl1);

	// 
	cv::Mat solidSkeletonImg;
	this->GetSolidSkeletonImg(solidSkeletonImg);
	cv::subtract(chnl0, solidSkeletonImg, chnl0);
	cv::subtract(chnl2, solidSkeletonImg, chnl2);
	cv::add(chnl1, solidSkeletonImg, chnl1);

	//
	cv::Mat voidSkeletonImg;
	this->GetVoidSkeletonImg(voidSkeletonImg);
	cv::subtract(chnl0, voidSkeletonImg, chnl0);
	cv::subtract(chnl1, voidSkeletonImg, chnl1);
	cv::add(chnl2, voidSkeletonImg, chnl2);

	channels[2] = chnl2.clone(); // r
	channels[1] = chnl1.clone(); // g
	channels[0] = chnl0.clone(); // b
	cv::merge(channels, background_img_);
}

void CImgData::convert_skeleton_to_img(bool is_solid, cv::Mat & img)
{
	int w, h;
	GetBoundingDomain(w, h);
	cv::Mat skel_img = cv::Mat(h, w, CV_8UC1);
	skel_img.setTo(0);
	
	if (is_solid) {
		for (int i = 0; i < solid_skeletal_pts_.size(); i++)
		{
			cv::Point p = solid_skeletal_pts_[i];
			skel_img.at<uchar>(p.y, p.x) = 255;
		}
	}
	else {
		for (int i = 0; i < void_skeletal_pts_.size(); i++)
		{
			cv::Point p = void_skeletal_pts_[i];
			skel_img.at<uchar>(p.y, p.x) = 255;
		}
	}

	skel_img.copyTo(img);
}

void CImgData::find_contours(std::vector<std::vector<double>>& edge_pts, const cv::Mat bw)
{
	//cv::Mat tmpBW;
	//bitwise_not(bw, bw);

	// get contours
	std::vector<std::vector<cv::Point> > contours;
	findContours(bw, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_NONE);

	for (int i = 0; i < contours.size(); ++i)
	{
		for (int j = 0; j < contours[i].size(); ++j)
		{
			std::vector<double> pt;
			pt.push_back(contours[i][j].x);
			pt.push_back(contours[i][j].y);
			edge_pts.push_back(pt);
		}
	}
}
