#include "cv_window.h"

CCVWindow::CCVWindow(char* fname):
	winName_("DRAW"), timeElapse_(0)
{
	src_img_ = cv::imread(fname);
	draw_ = cv::Mat(src_img_.size(), CV_8U);
}

CCVWindow::~CCVWindow()
{
}

void CCVWindow::setWindowSize(int _w, int _h)
{
	width_ = _w;
	height_ = _h;
}

void CCVWindow::setWindowParameters(char * winName, int timeElapse)
{
	winName_ = winName;
	timeElapse_ = timeElapse;
}

void CCVWindow::showVoronoiDiagram(CVoronoiDiagram * VD)
{
	VD->getVoronoiSegments(vsegments_);
	crop_voronoi_segments();
	visualize_diagram();
}

void CCVWindow::getMouseClickForAddingVSite(int &x, int &y)
{
	cv::Point *p;
	cv::setMouseCallback(winName_, onMouse, &p);
	// update VD
	WPoint wp(p->x, p->y, 0.0);
	VD_->addSite(wp);
}

void CCVWindow::crop_voronoi_segments()
{
	for (int i = 0; i < vsegments_.size(); i++)
	{
		std::pair<WPoint, WPoint> seg = vsegments_[i];
		// crop with [0,0->0,1]

		// crop with [0,1->1,1]

		// crop with [1,1->1,0]

		// crop with [1,0->0,0]

	}

	// finally we will have a updated list of voronoi segments for visualization
}

void CCVWindow::visualize_diagram()
{
	for (int i = 0; i < vsegments_.size(); i++)
	{
		std::pair<WPoint, WPoint> seg = vsegments_[i];
		// draw each segment
	}

	cv::namedWindow(winName_, cv::WINDOW_AUTOSIZE);
	cv::imshow(winName_, draw_);
	cv::waitKey(timeElapse_);
}

void CCVWindow::onMouse(int event, int x, int y, int flags, void * ptr)
{
	cv::Point*p = (cv::Point*)ptr;

	// when click left button, we record the position of the cursor for adding a new site
	if (event == cv::EVENT_LBUTTONDOWN)
	{
		std::cout << "Left button of the mouse is clicked - position (" << x << ", " << y << ")" << std::endl;
		p->x = x;
		p->y = y;
	}
	// when click right button, we record the position of the cursor for selecting a new site
	else if (event == cv::EVENT_RBUTTONDOWN)
	{
		std::cout << "Right button of the mouse is clicked - position (" << x << ", " << y << ")" << std::endl;
	}
	else if (event == cv::EVENT_MBUTTONDOWN)
	{
		std::cout << "Middle button of the mouse is clicked - position (" << x << ", " << y << ")" << std::endl;
	}
	else if (event == cv::EVENT_MOUSEMOVE)
	{
		std::cout << "Mouse move over the window - position (" << x << ", " << y << ")" << std::endl;
	}
}