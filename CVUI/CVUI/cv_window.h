#ifndef C_CV_WINDOW_H
#define C_CV_WINDOW_H

#include <iostream>

// INCLUDE OPENCV DIR
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include <opencv2/opencv.hpp>

#include "voronoi_drawer.h"
#include "../DataColle/voronoi_diagram.h"

struct cv_color_palette {
	cv::Scalar color_0 = cv::Scalar(13, 177, 240);
	cv::Scalar color_1 = cv::Scalar(115, 227, 9);
	cv::Scalar color_2 = cv::Scalar(255, 0, 193);
	cv::Scalar color_3 = cv::Scalar(249, 105, 199);
	cv::Scalar color_4 = cv::Scalar(0, 4, 255);
	cv::Scalar color_5 = cv::Scalar(184, 44, 44);
	cv::Scalar color_6 = cv::Scalar(244, 159, 89);
	cv::Scalar color_7 = cv::Scalar(255, 250, 139);
	cv::Scalar color_8 = cv::Scalar(56, 107, 48);
	cv::Scalar color_9 = cv::Scalar(236, 45, 74);
	cv::Scalar color_default = cv::Scalar(255, 0, 0);
};

class CCVWindow
{
public:
	CCVWindow();
	~CCVWindow();

	// set data
	void SetVoronoiDrawer(CVoronoiDrawer* p_VoroDrawer);

	// set view size
	void SetWindowSize(int _w, int _h);

	// set parameters
	void SetWindowParameters(char* winName, int timeElapse);

	// show the voronoi diagram (actually its segments) in the window
	void ShowUpdatedVoronoiDiagram(int elapseTime);

	// show the reconstruction points in terms of their clusters
	void ShowReconstructionPointClusters(std::vector<Point_2> rp_list, std::vector<int> rp_label_list, int elapseTime);
	
	//// detect mouse click that indicates adding voronoi site
	//// return a flag to indicate the mouse event
	//int getMouseClick(int &x, int &y);

private:
	//// source image
	//cv::Mat src_img_;

	//// source image with overlayed voronoi segments
	//cv::Mat draw_;

	// view size
	int width_;
	int height_;

	// parameter
	char* winName_;
	int timeElapse_;

	// VD drawer
	//CVoronoiDrawer voro_drawer_;
	CVoronoiDrawer* p_voro_drawer_;

	void drawSkeleton(int elapseTime, bool is_overlay = true);

	//// static functions
	//static void CCVWindow::on_mouse(int event, int x, int y, int flags, void * userdata);
	//static void CCVWindow::on_trackbar(int, void *);

	void drawLine(cv::Mat & img, cv::Point p1, cv::Point p2, int label, int lw = 1);
	void drawPolyLine(cv::Mat & img, std::vector<cv::Point> plist, int label, int lw = 1);
	cv::Point convert_to_cvPoint(Point_2);
};





#endif // !C_CV_WINDOW_H



