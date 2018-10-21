#include "cv_window.h"

CCVWindow::CCVWindow(char* fname):
	winName_("DRAW"), timeElapse_(0), 
	trackBarName_("Weight"), slider_max_(100), slider_(0)
{
	src_img_ = cv::imread(fname);
	draw_ = cv::Mat(src_img_.size(), CV_8U);
	p_VoroDrawer_ = new CVoronoiDrawer;
	p_cursorparameters_ = new CCursorParameters;

	cv::namedWindow(winName_, cv::WINDOW_AUTOSIZE);
	cv::createTrackbar(trackBarName_, winName_, &slider_, slider_max_, on_trackbar);

}

CCVWindow::~CCVWindow()
{
	delete p_cursorparameters_;
	delete p_VoroDrawer_;
	delete winName_;
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

void CCVWindow::showUpdatedVoronoiDiagram()
{
	std::vector<std::pair<cv::Point, cv::Point>> voronoi_edges;
	voronoi_edges = p_VoroDrawer_->DrawVoronoi();
	
	// draw voronoi_edges to the window
}

int CCVWindow::getMouseClick(int &x, int &y)
{
	cv::setMouseCallback(winName_, on_mouse, &p_cursorparameters_);
	cv::Point p = p_cursorparameters_->p_;
	int flag = p_cursorparameters_->flag_;

	//// left-click to add a new site to VD
	//if (flag == 1)
	//{
	//	// add point and update VD
	//	cv::Point p2(p.x, p.y);
	//	p_VD_->addSite(p2);
	//}
	//// right-click to add a new site to VD
	//else if (flag == 2)
	//{
	//	// do we need to reset the trackbar?
	//	cv::Point query(p.x, p.y);
	//	WPoint p2 = p_VoroDrawer_->pickSite(query);
	//	if (false)
	//	{
	//		// get weight from the trackbar
	//		double w = p2.weight();
	//		cv::setTrackbarPos(trackBarName_, winName_, (int)w*10.0);
	//		// can we have some update?
	//		slider_ = cv::getTrackbarPos(trackBarName_, winName_);
	//		w = (double)slider_ / 10.0;
	//		p_VoroDrawer_->setWeightToPickedSite(p2.point(), w);
	//	}
	//}

	return flag;
}

void CCVWindow::on_mouse(int event, int x, int y, int flags, void * ptr)
{
	CCursorParameters* cp = (CCursorParameters*)ptr;
	// when click left button, we record the position of the cursor for adding a new site
	if (event == cv::EVENT_LBUTTONDOWN)
	{
		std::cout << "Left button of the mouse is clicked - position (" << x << ", " << y << ")" << std::endl;
		cp->p_.x = x;
		cp->p_.y = y;
		cp->flag_ = flags;
	}
	// when click right button, we record the position of the cursor for selecting a new site
	else if (event == cv::EVENT_RBUTTONDOWN)
	{
		std::cout << "Right button of the mouse is clicked - position (" << x << ", " << y << ")" << std::endl;
		cp->p_.x = x;
		cp->p_.y = y;
		cp->flag_ = flags;
	}
	//else if (event == cv::EVENT_MBUTTONDOWN)
	//{
	//	std::cout << "Middle button of the mouse is clicked - position (" << x << ", " << y << ")" << std::endl;
	//}
	//else if (event == cv::EVENT_MOUSEMOVE)
	//{
	//	std::cout << "Mouse move over the window - position (" << x << ", " << y << ")" << std::endl;
	//}
}

void CCVWindow::on_trackbar(int, void *)
{
	//alpha = (double)alpha_slider / alpha_slider_max;
	//beta = (1.0 - alpha);

	//addWeighted(src1, alpha, src2, beta, 0.0, dst);

	//imshow("Linear Blend", dst);
}
