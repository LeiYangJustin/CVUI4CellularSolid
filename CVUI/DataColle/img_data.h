#ifndef C_IMG_DATA

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
	CImgData();
	~CImgData();
	bool ReadImgFromFile(const char* filename);
	void SetSrcImg(cv::Mat src_img);

	cv::Mat GetVoidImg() { return void_img_; };
	cv::Mat GetSolidImg() { return solid_img_; };
	cv::Mat GetSrcImg() { return src_img_; };
	void getBoundingDomain(int & width, int & height);

private:
	cv::Mat src_img_;
	cv::Mat solid_img_;
	cv::Mat void_img_;	

	void binarize_src_img();
};

#endif // !C_IMG_DATA




