#include "cv_window.h"

struct tmp_Data
{
	tmp_Data(std::vector<std::pair<cv::Point, cv::Point>> cv_voronoi_edges,
		std::map<int, std::vector<cv::Point>> map_lid_to_rplist)
	{
		cv_voronoi_edges_ = cv_voronoi_edges;
		map_lid_to_rplist_ = map_lid_to_rplist;
	}
	std::vector<std::pair<cv::Point, cv::Point>> cv_voronoi_edges_;
	std::map<int, std::vector<cv::Point>> map_lid_to_rplist_;
	double compute_dist_from_point_to_segment(cv::Point p, cv::Point src, cv::Point tgt)
	{
		return compute_dist_from_point_to_segment(Point_2(p.x, p.y), Segment_2(Point_2(src.x, src.y), Point_2(tgt.x, tgt.y)));
	};
	double compute_dist_from_point_to_segment(Point_2 p, Segment_2 e)
	{
		Vector_2 a = p - e.source();
		Vector_2 v = e.target() - e.source();
		Vector_2 proj_a_v = v*(a*v / (v*v));
		double t = sqrt(proj_a_v.squared_length()) / sqrt(v.squared_length());
		
		double s = a*v;
		if (s < 0)
			t = -t;
		if (t >= 0 && t <= 1)
		{
			Vector_2 d = a - proj_a_v;
			return sqrt(d.squared_length());
		}
		else if (t < 0) {
			return sqrt(a.squared_length());
		}
		else {
			return sqrt((p - e.target()).squared_length());
		}
	};
};

static void onMouse(int event, int x, int y, int flag, void* ptr)
{
	if (event != cv::EVENT_LBUTTONDOWN)
		return;

	tmp_Data *p = (tmp_Data*)ptr;

	// cv::Point* p = (cv::Point*) ptr;
	//p->x = x;
	//p->y = y;
	//std::cout << "click at " << *p << std::endl;

	int min_id = -1;
	double min_dist = 1000;
	for (auto  miter = p->map_lid_to_rplist_.begin(); miter != p->map_lid_to_rplist_.end(); ++miter)
	{
		std::vector<cv::Point> plist = miter->second;
		for (int i = 0; i < plist.size(); i++)
		{
			double dist = (Point_2(x, y) - Point_2(plist[i].x, plist[i].y)).squared_length();
			if (dist < min_dist)
			{
				min_dist = dist;
				min_id = miter->first;
			}
		}
	}
	std::cout << "click at p" << cv::Point(x,y) << " with label " << min_id
		<< " whose corresponding segment is " 
		<< p->cv_voronoi_edges_[min_id].first
		<< " -> " 
		<< p->cv_voronoi_edges_[min_id].second
		<< std::endl;
	std::cout << "The dist from the click to the segment is " <<
		p->compute_dist_from_point_to_segment(cv::Point(x, y),
			p->cv_voronoi_edges_[min_id].first,
			p->cv_voronoi_edges_[min_id].second)
		<< std::endl;
	std::cout << std::endl;
		
}

CCVWindow::CCVWindow():
	winName_("DRAW"), timeElapse_(0)
{
	//p_VoroDrawer_ = new CVoronoiDrawer;
	p_voro_drawer_ = new CVoronoiDrawer;

}

CCVWindow::~CCVWindow()
{
	//delete p_cursorparameters_;
	delete p_voro_drawer_;
}

void CCVWindow::SetVoronoiDrawer(CVoronoiDrawer* p_voro_drawer)
{
	p_voro_drawer_ = p_voro_drawer;
}

void CCVWindow::SetWindowSize(int _w, int _h)
{
	width_ = _w;
	height_ = _h;
}

void CCVWindow::SetWindowParameters(char * winName, int timeElapse)
{
	winName_ = winName;
	timeElapse_ = timeElapse;
}

void CCVWindow::ShowUpdatedVoronoiDiagram(int elapseTime)
{
	drawSkeleton(elapseTime*1000 /*make it miliseconds*/);
}

void CCVWindow::ShowReconstructionPointClusters(std::vector<Point_2> rp_list, std::vector<int> rp_label_list, int elapseTime)
{
	// draw the background
	CImgData* p_img_data = p_voro_drawer_->GetImgData();

	cv::Mat result_img;
	std::vector<cv::Mat> channels(3);
	cv::Mat chnl2 = p_img_data->GetSolidImg().clone();
	cv::Mat chnl1 = p_img_data->GetSolidImg().clone();
	cv::Mat chnl0 = p_img_data->GetSolidImg().clone();
	channels[2] = chnl2.clone(); // r
	channels[1] = chnl1.clone(); // g
	channels[0] = chnl0.clone(); // b
	cv::merge(channels, result_img);

	std::map<int, std::vector<cv::Point>> map_lid_to_rplist;
	for (int i = 0; i != rp_list.size(); ++i)
	{
		// add the point to the point set
		if (map_lid_to_rplist.find(rp_label_list[i]) != map_lid_to_rplist.end())
		{
			map_lid_to_rplist[rp_label_list[i]].push_back(convert_to_cvPoint(rp_list[i]));
		}
		else {
			std::vector<cv::Point> emptyList;
			emptyList.push_back(convert_to_cvPoint(rp_list[i]));
			map_lid_to_rplist.insert(std::make_pair(rp_label_list[i], emptyList));
		}
	}

	std::vector<std::pair<cv::Point, cv::Point>> cv_voronoi_edges = p_voro_drawer_->GetVoronoiEdges();

	for (auto miter = map_lid_to_rplist.begin(); miter != map_lid_to_rplist.end(); ++miter)
	{
		if (miter->first % 1 == 0)
		{
			drawPolyLine(result_img, miter->second, miter->first, 2);
			drawLine(result_img,
				cv_voronoi_edges[miter->first].first, cv_voronoi_edges[miter->first].second,
				miter->first);
		}
	}

	//for (int i = 0; i < rp_list.size(); i++)
	//{
	//	drawLine(result_img,
	//		cv::Point(rp_list[i].x(), rp_list[i].y()), cv::Point(rp_list[i].x(), rp_list[i].y()),
	//		rp_label_list[i]);
	//}
	
	cv::namedWindow("Reconstruction", cv::WINDOW_AUTOSIZE);
	cv::imshow("Reconstruction", result_img);

	tmp_Data p(cv_voronoi_edges, map_lid_to_rplist);
	cv::setMouseCallback("Reconstruction", onMouse, (void*)&p);
	std::cout << "Please Pick a point!";

	for (;;)
	{
		int c = cv::waitKey(0);
		if (c == 27)
		{
			std::cout << "Exit ...\n";
			break;
		}
	}
	cv::destroyAllWindows();
}

void CCVWindow::drawSkeleton(int elapseTime, bool is_overlay)
{
	// draw the background
	CImgData* p_img_data = p_voro_drawer_->GetImgData();
	// OVERLAY IMAGE
	cv::Mat overlayImg(p_img_data->GetSolidImg().size(), CV_8UC3, cv::Scalar(0, 0, 0));
	if (is_overlay) 
	{
		std::vector<cv::Mat> channels(3);
		cv::Mat chnl2 = p_img_data->GetSolidImg().clone();
		cv::Mat chnl1 = p_img_data->GetSolidImg().clone();
		cv::Mat chnl0 = p_img_data->GetSolidImg().clone();

		////
		//cv::subtract(chnl0, skelImg_secondary_, chnl0);
		//cv::subtract(chnl2, skelImg_secondary_, chnl2);
		//cv::add(chnl1, skelImg_secondary_, chnl1);

		// 
		cv::Mat solidSkeletonImg;
		p_img_data->GetSolidSkeletonImg(solidSkeletonImg);
		cv::subtract(chnl0, solidSkeletonImg, chnl0);
		cv::subtract(chnl2, solidSkeletonImg, chnl2);
		cv::add(chnl1, solidSkeletonImg, chnl1);

		//
		cv::Mat voidSkeletonImg;
		p_img_data->GetVoidSkeletonImg(voidSkeletonImg);
		cv::subtract(chnl0, voidSkeletonImg, chnl0);
		cv::subtract(chnl1, voidSkeletonImg, chnl1);
		cv::add(chnl2, voidSkeletonImg, chnl2);

		channels[2] = chnl2.clone(); // r
		channels[1] = chnl1.clone(); // g
		channels[0] = chnl0.clone(); // b
		cv::merge(channels, overlayImg);

		cv::imshow("material with two skeletons", overlayImg);
		cv::waitKey(5);
	}

	cv::Mat result_img;
	std::vector<cv::Mat> channels(3);
	cv::Mat chnl2 = p_img_data->GetSolidImg().clone();
	cv::Mat chnl1 = p_img_data->GetSolidImg().clone();
	cv::Mat chnl0 = p_img_data->GetSolidImg().clone();
	channels[2] = chnl2.clone(); // r
	channels[1] = chnl1.clone(); // g
	channels[0] = chnl0.clone(); // b
	cv::merge(channels, result_img);

	// draw voronoi edges
	std::vector<std::pair<cv::Point, cv::Point>> voronoi_edges;
	voronoi_edges = p_voro_drawer_->GetVoronoiEdges();
	for (int i = 0; i < voronoi_edges.size(); i++) {
		//cv::line(overlayImg, voronoi_edges[i].first, voronoi_edges[i].second, cv::Scalar(255, 0, 0), 2, 2);
		//drawLine(overlayImg, voronoi_edges[i].first, voronoi_edges[i].second, -1, 1);
	}

	// draw fitted edges
	std::vector<std::pair<cv::Point, cv::Point>> fitted_edges;
	fitted_edges = p_voro_drawer_->GetFittedEdges();
	for (int i = 0; i < fitted_edges.size(); i++) {
		//cv::line(overlayImg, fitted_edges[i].first, fitted_edges[i].second, cv::Scalar(0, 255, 0), 2, 2);
		drawLine(overlayImg, fitted_edges[i].first, fitted_edges[i].second, -1, 2);
	}

	//// draw fitting base pts
	//std::vector<std::pair<cv::Point, cv::Point>> fitting_base_pts;
	//fitting_base_pts = p_voro_drawer_->GetFittingBasePoints();
	//for (int i = 0; i < fitting_base_pts.size(); i++) {
	//	//cv::line(result_img, fitting_base_pts[i].first, fitting_base_pts[i].second, cv::Scalar(0, 0, 255), 2, 2);
	//	drawLine(overlayImg, fitting_base_pts[i].first, fitting_base_pts[i].second, i, 2);
	//}
	// draw fitting base pts
	std::vector<std::vector<cv::Point>> fitting_base_pts;
	fitting_base_pts = p_voro_drawer_->GetFittingBasePointsList();
	for (int i = 0; i < fitting_base_pts.size(); i++) {
		//cv::line(result_img, fitting_base_pts[i].first, fitting_base_pts[i].second, cv::Scalar(0, 0, 255), 2, 2);
		drawPolyLine(overlayImg, fitting_base_pts[i], i, 2);
	}

	//cv::namedWindow(winName_, cv::WINDOW_KEEPRATIO);
	cv::imshow(winName_, overlayImg);
	cv::waitKey(elapseTime);

	//if (b_output_)
	//{
	//	outputImgFile("result", char_fname_, windowName, overlayImg);
	//	//char fn_output[128];
	//	//std::strcpy(fn_output, "result/");
	//	//std::strcat(fn_output, char_fname_);
	//	//std::strcat(fn_output, "/");
	//	//std::strcat(fn_output, "skeleton");
	//	//std::strcat(fn_output, char_ext_);
	//	//printf("%s\n", fn_output);
	//	//cv::imwrite(fn_output, overlayImg);
	//}
}

void CCVWindow::drawLine(cv::Mat & img, cv::Point p1, cv::Point p2, int label, int lw)
{
	cv_color_palette color_pat;
	switch (label % 10)
	{
	case 0:
		cv::arrowedLine(img, p1, p2, color_pat.color_0, lw, lw);
		break;
	case 1:
		cv::arrowedLine(img, p1, p2, color_pat.color_1, lw, lw);
		break;
	case 2:
		cv::arrowedLine(img, p1, p2, color_pat.color_2, lw, lw);
		break;
	case 3:
		cv::arrowedLine(img, p1, p2, color_pat.color_3, lw, lw);
		break;
	case 4:
		cv::arrowedLine(img, p1, p2, color_pat.color_4, lw, lw);
		break;
	case 5:
		cv::arrowedLine(img, p1, p2, color_pat.color_5, lw, lw);
		break;
	case 6:
		cv::arrowedLine(img, p1, p2, color_pat.color_6, lw, lw);
		break;
	case 7:
		cv::arrowedLine(img, p1, p2, color_pat.color_7, lw, lw);
		break;
	case 8:
		cv::arrowedLine(img, p1, p2, color_pat.color_8, lw, lw);
		break;
	case 9:
		cv::arrowedLine(img, p1, p2, color_pat.color_9, lw, lw);
		break;
	default:
		cv::arrowedLine(img, p1, p2, color_pat.color_default, lw, lw);
		break;
	}
}

void CCVWindow::drawPolyLine(cv::Mat & img, std::vector<cv::Point> plist, int label, int lw)
{
	cv_color_palette color_pat;
	switch (label % 10)
	{
	case 0:
		cv::polylines(img, plist, false, color_pat.color_0, lw);
		break;
	case 1:
		cv::polylines(img, plist, false, color_pat.color_1, lw);
		break;
	case 2:
		cv::polylines(img, plist, false, color_pat.color_2, lw);
		break;
	case 3:
		cv::polylines(img, plist, false, color_pat.color_3, lw);
		break;
	case 4:
		cv::polylines(img, plist, false, color_pat.color_4, lw);
		break;
	case 5:
		cv::polylines(img, plist, false, color_pat.color_5, lw);
		break;
	case 6:
		cv::polylines(img, plist, false, color_pat.color_6, lw);
		break;
	case 7:
		cv::polylines(img, plist, false, color_pat.color_7, lw);
		break;
	case 8:
		cv::polylines(img, plist, false, color_pat.color_8, lw);
		break;
	case 9:
		cv::polylines(img, plist, false, color_pat.color_9, lw);
		break;
	default:
		cv::polylines(img, plist, false, color_pat.color_default, lw);
		break;
	}

	
}

cv::Point CCVWindow::convert_to_cvPoint(Point_2 p)
{
	return cv::Point(p.x(), p.y());
}

//
//int CCVWindow::getMouseClick(int &x, int &y)
//{
//	cv::setMouseCallback(winName_, on_mouse, &p_cursorparameters_);
//	cv::Point p = p_cursorparameters_->p_;
//	int flag = p_cursorparameters_->flag_;
//
//	//// left-click to add a new site to VD
//	//if (flag == 1)
//	//{
//	//	// add point and update VD
//	//	cv::Point p2(p.x, p.y);
//	//	p_VD_->addSite(p2);
//	//}
//	//// right-click to add a new site to VD
//	//else if (flag == 2)
//	//{
//	//	// do we need to reset the trackbar?
//	//	cv::Point query(p.x, p.y);
//	//	WPoint p2 = p_VoroDrawer_->pickSite(query);
//	//	if (false)
//	//	{
//	//		// get weight from the trackbar
//	//		double w = p2.weight();
//	//		cv::setTrackbarPos(trackBarName_, winName_, (int)w*10.0);
//	//		// can we have some update?
//	//		slider_ = cv::getTrackbarPos(trackBarName_, winName_);
//	//		w = (double)slider_ / 10.0;
//	//		p_VoroDrawer_->setWeightToPickedSite(p2.point(), w);
//	//	}
//	//}
//
//	return flag;
//}
//
//void CCVWindow::on_mouse(int event, int x, int y, int flags, void * ptr)
//{
//	CCursorParameters* cp = (CCursorParameters*)ptr;
//	// when click left button, we record the position of the cursor for adding a new site
//	if (event == cv::EVENT_LBUTTONDOWN)
//	{
//		std::cout << "Left button of the mouse is clicked - position (" << x << ", " << y << ")" << std::endl;
//		cp->p_.x = x;
//		cp->p_.y = y;
//		cp->flag_ = flags;
//	}
//	// when click right button, we record the position of the cursor for selecting a new site
//	else if (event == cv::EVENT_RBUTTONDOWN)
//	{
//		std::cout << "Right button of the mouse is clicked - position (" << x << ", " << y << ")" << std::endl;
//		cp->p_.x = x;
//		cp->p_.y = y;
//		cp->flag_ = flags;
//	}
//	//else if (event == cv::EVENT_MBUTTONDOWN)
//	//{
//	//	std::cout << "Middle button of the mouse is clicked - position (" << x << ", " << y << ")" << std::endl;
//	//}
//	//else if (event == cv::EVENT_MOUSEMOVE)
//	//{
//	//	std::cout << "Mouse move over the window - position (" << x << ", " << y << ")" << std::endl;
//	//}
//}
//
//void CCVWindow::on_trackbar(int, void *)
//{
//	//alpha = (double)alpha_slider / alpha_slider_max;
//	//beta = (1.0 - alpha);
//
//	//addWeighted(src1, alpha, src2, beta, 0.0, dst);
//
//	//imshow("Linear Blend", dst);
//}
