#ifndef C_CV_WINDOW

#include "../DataColle/voronoi_diagram.h"
#include <iostream>
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"

class CCursorParameters
{
public:
	cv::Point p_;
	int flag_;
};

class CCVWindow
{
public:
	CCVWindow(char* fname);
	~CCVWindow();

	// set view size
	void setWindowSize(int _w, int _h);

	// set parameters
	void setWindowParameters(char* winName, int timeElapse);

	// show the voronoi diagram (actually its segments) in the window
	void showUpdatedVoronoiDiagram();

	// detect mouse click that indicates adding voronoi site
	// return a flag to indicate the mouse event
	int getMouseClick(int &x, int &y);

private:
	void crop_voronoi_segments();
	void visualize_diagram();

	// static functions
	static void CCVWindow::on_mouse(int event, int x, int y, int flags, void * userdata);
	static void CCVWindow::on_trackbar(int, void *);

	// private variables
	std::vector<std::pair<iPoint2, iPoint2>> vsegments_;

	// source image
	cv::Mat src_img_;

	// source image with overlayed voronoi segments
	cv::Mat draw_;

	// view size
	int width_;
	int height_;

	// parameter
	char* winName_;
	int timeElapse_;
	char* trackBarName_;
	int slider_;
	int slider_max_;
	// VD
	CVoronoiDiagram* p_VD_;

	// cursor parameter
	CCursorParameters* p_cursorparameters_;
};

#endif // !C_CV_WINDOW



