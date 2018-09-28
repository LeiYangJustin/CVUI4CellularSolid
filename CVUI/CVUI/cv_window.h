#ifndef C_CV_WINDOW

#include "../DataColle/voronoi_diagram.h"
#include <iostream>
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"

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
	void showVoronoiDiagram(CVoronoiDiagram *VD);

	// detect mouse click that indicates adding voronoi site
	void getMouseClickForAddingVSite(int &x, int &y);

private:
	void crop_voronoi_segments();
	void visualize_diagram();
	static void CCVWindow::onMouse(int event, int x, int y, int flags, void * userdata);

	// private variables
	std::vector<std::pair<WPoint, WPoint>> vsegments_;

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

	// VD
	CVoronoiDiagram* VD_;
};

#endif // !C_CV_WINDOW



