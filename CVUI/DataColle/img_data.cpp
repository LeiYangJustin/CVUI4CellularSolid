#include "img_data.h"



CImgData::CImgData()
{
}


CImgData::~CImgData()
{
}

void CImgData::SetSrcImg(cv::Mat src_img)
{
#ifdef _DEBUG
	int dRows = 101;
	//int dCols = 101;
	std::cout << "We are in the debugging mode" << std::endl;
#else 
	int dRows = 401;
	//int dCols = 501;
#endif
	// resize the image
	int dCols = double(dRows) / double(src_img.rows)*double(src_img.cols);
	cv::resize(src_img, src_img_, cv::Size(dCols, dRows));

	// BINARIZE THE EXEMPLAR
	binarize_src_img();
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


void CImgData::getBoundingDomain(int & width, int & height)
{
	width = src_img_.cols;
	height = src_img_.rows;
}

void CImgData::binarize_src_img()
{
	// BINARIZE THE EXEMPLAR
	cv::Mat grayExplImg(src_img_.size(), CV_8UC1);
	cv::cvtColor(src_img_, grayExplImg, cv::COLOR_BGR2GRAY);
	cv::blur(grayExplImg, solid_img_, cv::Size(3, 3));
	cv::threshold(solid_img_, solid_img_, 50, 255, cv::THRESH_BINARY | cv::THRESH_OTSU);
	cv::bitwise_not(solid_img_, void_img_);
}
