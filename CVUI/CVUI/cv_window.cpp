#include "cv_window.h"

static double compute_dist_from_point_to_segment(Point_2 p, Segment_2 e)
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

static double compute_dist_from_point_to_segment(cv::Point p, cv::Point src, cv::Point tgt)
{
	return compute_dist_from_point_to_segment(Point_2(p.x, p.y), Segment_2(Point_2(src.x, src.y), Point_2(tgt.x, tgt.y)));
};


static double two_cvpoint_distance(cv::Point p1, cv::Point p2)
{
	cv::Point d = p1 - p2;
	return sqrt(d.dot(d));
}

struct UI_DATA_DRAW_SKELETON
{
	std::string winname_;
	std::vector<std::pair<cv::Point, cv::Point>> cv_voronoi_edges_;
	std::vector<std::pair<cv::Point, cv::Point>> cv_fitted_edges_;
	CVoronoiDiagram* pVD_;

	UI_DATA_DRAW_SKELETON(
		std::vector<std::pair<cv::Point, cv::Point>> cv_voronoi_edges,
		std::vector<std::pair<cv::Point, cv::Point>> cv_fitted_edges,
		CVoronoiDiagram* pVD,
		std::string winname)
	{
		cv_voronoi_edges_ = cv_voronoi_edges;
		cv_fitted_edges_ = cv_fitted_edges;
		pVD_ = pVD;
		winname_ = winname;
	}
};

struct tmp_Data
{
	std::string winname_;
	std::vector<std::pair<cv::Point, cv::Point>> cv_voronoi_edges_;
	std::map<int, std::vector<cv::Point>> map_lid_to_rplist_;

	tmp_Data(std::vector<std::pair<cv::Point, cv::Point>> cv_voronoi_edges,
		std::map<int, std::vector<cv::Point>> map_lid_to_rplist, 
		std::string winname)
	{
		cv_voronoi_edges_ = cv_voronoi_edges;
		map_lid_to_rplist_ = map_lid_to_rplist;
		winname_ = winname;
	}
};

static void onMouseSkeleton(int event, int x, int y, int flag, void* ptr)
{
	if (event != cv::EVENT_LBUTTONDOWN)
		return;

	UI_DATA_DRAW_SKELETON *p = (UI_DATA_DRAW_SKELETON*)ptr;

	int min_id_vedges = -1, id = 0;
	double min_dist_vedges = 1000;
	for (auto viter = p->cv_voronoi_edges_.begin(); viter != p->cv_voronoi_edges_.end(); ++viter, ++id)
	{
		double dist = compute_dist_from_point_to_segment(cv::Point(x, y), viter->first, viter->second);
		if (dist < min_dist_vedges) {
			min_dist_vedges = dist;
			min_id_vedges = id;
		}
	}

	int min_id_dedges = -1;
	id = 0;
	double min_dist_dedges = 1000;
	for (auto viter = p->cv_fitted_edges_.begin(); viter != p->cv_fitted_edges_.end(); ++viter, ++id)
	{
		double dist = compute_dist_from_point_to_segment(cv::Point(x, y), viter->first, viter->second);
		if (dist < min_dist_dedges) {
			min_dist_dedges = dist;
			min_id_dedges = id;
		}
	}

	std::cout << "Window Name: " << p->winname_ << std::endl;
	std::cout << "click at p" << cv::Point(x, y)
		<< " whose corresponding segment is "
		<< p->cv_voronoi_edges_[min_id_vedges].first
		<< " -> "
		<< p->cv_voronoi_edges_[min_id_vedges].second
		<< std::endl;
	std::cout << " and whose corresponding segment is "
		<< p->cv_fitted_edges_[min_id_dedges].first
		<< " -> "
		<< p->cv_fitted_edges_[min_id_dedges].second
		<< std::endl;

	//typedef Regular_triangulation::Finite_edges_iterator RT_Finite_Edge_Iter;
	//typedef Regular_triangulation::Face_handle FaceHandle;

	//Regular_triangulation* rt = p->pVD_->getTriangulation();
	//FaceHandle fh1 = rt->locate(WPoint(p->cv_voronoi_edges_[min_id_vedges].first.x, p->cv_voronoi_edges_[min_id_vedges].first.y));
	//FaceHandle fh2 = rt->locate(WPoint(p->cv_voronoi_edges_[min_id_vedges].second.x, p->cv_voronoi_edges_[min_id_vedges].second.y));

	//FaceHandle src_fh = NULL, tgt_fh = NULL;
	//double src_dist = 1000, tgt_dist = 1000;
	//for (auto fiter = rt->finite_faces_begin(); fiter != rt->finite_faces_end(); ++fiter)
	//{
	//	cv::Point cvp = cv::Point(rt->dual(fiter).x(), rt->dual(fiter).y());
	//	if (two_cvpoint_distance(cvp, p->cv_voronoi_edges_[min_id_vedges].first) < src_dist)
	//	{
	//		src_dist = two_cvpoint_distance(cvp, p->cv_voronoi_edges_[min_id_vedges].first);
	//		src_fh = fiter;
	//	}
	//	if (two_cvpoint_distance(cvp, p->cv_voronoi_edges_[min_id_vedges].second) < tgt_dist)
	//	{
	//		tgt_dist = two_cvpoint_distance(cvp, p->cv_voronoi_edges_[min_id_vedges].second);
	//		tgt_fh = fiter;
	//	}
	//}
	//if (src_fh != NULL && tgt_fh != NULL)
	//{
	//	std::cout << "The dual edge is: "
	//		<< cv::Point(rt->dual(src_fh).x(), rt->dual(src_fh).y())
	//		<< " -> "
	//		<< cv::Point(rt->dual(tgt_fh).x(), rt->dual(tgt_fh).y())
	//		<< std::endl;

	//	Vector_2 src_pt(0, 0);
	//	Vector_2 tgt_pt(0, 0);
	//	std::vector<Point_2> src_plist = src_fh->info();
	//	for (int i = 0; i < src_plist.size(); i++)
	//	{
	//		src_pt += (src_plist[i] - Point_2(0, 0));
	//	}
	//	src_pt /= double(src_plist.size());

	//	std::vector<Point_2> tgt_plist = tgt_fh->info();
	//	for (int i = 0; i < tgt_plist.size(); i++)
	//	{
	//		tgt_pt += (tgt_plist[i] - Point_2(0, 0));
	//	}
	//	tgt_pt /= double(tgt_plist.size());

	//	std::cout << " whose fitted segment as computed is "
	//		<< cv::Point(src_pt.x(), src_pt.y())
	//		<< " -> "
	//		<< cv::Point(tgt_pt.x(), tgt_pt.y())
	//		<< std::endl;

	//	std::cout << "They, respectively, are based on:" << std::endl;
	//	std::copy(src_fh->info().begin(), src_fh->info().end(), std::ostream_iterator<Point_2>(std::cout, " / "));
	//	std::cout << std::endl;
	//	std::copy(tgt_fh->info().begin(), tgt_fh->info().end(), std::ostream_iterator<Point_2>(std::cout, " / "));
	//	std::cout << std::endl;
	//}

	//std::cout << " whose fitted segment is "
	//	<< cv::Point(fh_average_point.x(), fh_average_point.y())
	//	<< " -> "
	//	<< cv::Point(fh_oppo_average_point.x(), fh_oppo_average_point.y())
	//	<< std::endl;
	//std::cout << "The dist from the click to the segment is " <<
	//	compute_dist_from_point_to_segment(cv::Point(x, y),
	//		p->cv_voronoi_edges_[min_id].first,
	//		p->cv_voronoi_edges_[min_id].second)
	//	<< std::endl;
	//std::cout << std::endl;
}


static void onMouseReconstruction(int event, int x, int y, int flag, void* ptr)
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
	std::cout << "Window Name: " << p->winname_ << std::endl;
	std::cout << "click at p" << cv::Point(x,y) << " with label " << min_id
		<< " whose corresponding segment is " 
		<< p->cv_voronoi_edges_[min_id].first
		<< " -> " 
		<< p->cv_voronoi_edges_[min_id].second
		<< std::endl;
	std::cout << "The dist from the click to the segment is " <<
		compute_dist_from_point_to_segment(cv::Point(x, y),
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

	//for (auto miter = map_lid_to_rplist.begin(); miter != map_lid_to_rplist.end(); ++miter)
	//{
	//	if (miter->first == -1)
	//		continue;
	//	if (miter->first % 1 == 0)
	//	{
	//		drawPolyLine(result_img, miter->second, miter->first, 1);
	//		drawLine(result_img,
	//			cv_voronoi_edges[miter->first].first, cv_voronoi_edges[miter->first].second,
	//			miter->first, 2);
	//	}
	//}

	//for (int i = 0; i < cv_voronoi_edges.size(); i++)
	//{
	//	drawLine(result_img,
	//		cv_voronoi_edges[i].first, cv_voronoi_edges[i].second,
	//		i, 2);
	//}

	for (int i = 0; i < rp_list.size(); i++)
	{
		drawLine(result_img,
			cv::Point(rp_list[i].x(), rp_list[i].y()), cv::Point(rp_list[i].x(), rp_list[i].y()),
			rp_label_list[i], 2);
	}
	
	std::string winname = "Reconstruction";
	cv::namedWindow(winname, cv::WINDOW_AUTOSIZE);
	cv::imshow(winname, result_img);
	cv::waitKey(3);
	//tmp_Data p(cv_voronoi_edges, map_lid_to_rplist, winname);
	//cv::setMouseCallback(winname, onMouseReconstruction, (void*)&p);
	//std::cout << "Please Pick a point!";

	//for (;;)
	//{
	//	int c = cv::waitKey(0);
	//	if (c == 27)
	//	{
	//		std::cout << "Exit ...\n";
	//		break;
	//	}
	//}
	//cv::destroyAllWindows();
}

void CCVWindow::ShowUpdatedVoronoiDiagram(int elapseTime)
{
	// draw the background
	CImgData* p_img_data = p_voro_drawer_->GetImgData();
	// OVERLAY IMAGE
	cv::Mat overlayImg(p_img_data->GetSolidImg().size(), CV_8UC3, cv::Scalar(0, 0, 0));

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

	// draw voronoi edges
	std::vector<std::pair<cv::Point, cv::Point>> voronoi_edges;
	voronoi_edges = p_voro_drawer_->GetVoronoiEdges();
	for (int i = 0; i < voronoi_edges.size(); i++) {
		//cv::line(overlayImg, voronoi_edges[i].first, voronoi_edges[i].second, cv::Scalar(255, 0, 0), 2, 2);
		drawLine(overlayImg, voronoi_edges[i].first, voronoi_edges[i].second, -1, 1);
	}
	// draw voronoi edges
	std::vector<std::pair<cv::Point, cv::Point>> traingulate_edges;
	traingulate_edges = p_voro_drawer_->GetTriangulationEdges();
	for (int i = 0; i < traingulate_edges.size(); i++) {
		//cv::line(overlayImg, voronoi_edges[i].first, voronoi_edges[i].second, cv::Scalar(255, 0, 0), 2, 2);
		drawLine(overlayImg, traingulate_edges[i].first, traingulate_edges[i].second, 2, 1);
	}

	//// draw fitted edges
	//std::vector<std::pair<cv::Point, cv::Point>> fitted_edges;
	//fitted_edges = p_voro_drawer_->GetFittedEdges();
	////fitted_edges = p_voro_drawer_->GetDualEdges();
	//for (int i = 0; i < fitted_edges.size(); i++) {
	//	//cv::line(overlayImg, fitted_edges[i].first, fitted_edges[i].second, cv::Scalar(0, 255, 0), 4, 4);
	//	drawLine(overlayImg, fitted_edges[i].first, fitted_edges[i].second, -2, 2);
	//}

	cv::namedWindow(winName_);
	//cv::namedWindow(winName_, cv::WINDOW_KEEPRATIO);
	cv::imshow(winName_, overlayImg);
	cv::waitKey(elapseTime);

	//CVoronoiDiagram* pVD = p_voro_drawer_->GetVD();

	//UI_DATA_DRAW_SKELETON p(voronoi_edges, fitted_edges, pVD, winName_);
	//cv::setMouseCallback(winName_, onMouseSkeleton, (void*)&p);
	//std::cout << "Please Pick a point!";

	//for (;;)
	//{
	//	int c = cv::waitKey(0);
	//	if (c == 27)
	//	{
	//		std::cout << "Exit ...\n";
	//		break;
	//	}
	//}
	//cv::destroyAllWindows();


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

void CCVWindow::ShowUpdatedRegularTriangulation(int elapseTime)
{
	// draw the background
	CImgData* p_img_data = p_voro_drawer_->GetImgData();
	// OVERLAY IMAGE
	cv::Mat overlayImg(p_img_data->GetSolidImg().size(), CV_8UC3, cv::Scalar(0, 0, 0));

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

	// draw voronoi edges
	std::vector<std::pair<cv::Point, cv::Point>> edges;
	edges = p_voro_drawer_->GetTriangulationEdges();
	for (int i = 0; i < edges.size(); i++) {
		drawLine(overlayImg, edges[i].first, edges[i].second, -1, 1);
	}

	std::string winName = "Regular_Triangulation";
	cv::namedWindow(winName);
	cv::imshow(winName, overlayImg);
	cv::waitKey(elapseTime);
}

void CCVWindow::drawLine(cv::Mat & img, cv::Point p1, cv::Point p2, int label, int lw)
{
	cv_color_palette color_pat;
	float tipLength = 0.5;
	switch (label % 10)
	{
	case 0:
		cv::arrowedLine(img, p1, p2, color_pat.color_0, lw, lw, tipLength);
		break;
	case 1:
		cv::arrowedLine(img, p1, p2, color_pat.color_1, lw, lw, tipLength);
		break;
	case 2:
		cv::arrowedLine(img, p1, p2, color_pat.color_2, lw, lw, tipLength);
		break;
	case 3:
		cv::arrowedLine(img, p1, p2, color_pat.color_3, lw, lw, tipLength);
		break;
	case 4:
		cv::arrowedLine(img, p1, p2, color_pat.color_4, lw, lw, tipLength);
		break;
	case 5:
		cv::arrowedLine(img, p1, p2, color_pat.color_5, lw, lw, tipLength);
		break;
	case 6:
		cv::arrowedLine(img, p1, p2, color_pat.color_6, lw, lw, tipLength);
		break;
	case 7:
		cv::arrowedLine(img, p1, p2, color_pat.color_7, lw, lw, tipLength);
		break;
	case 8:
		cv::arrowedLine(img, p1, p2, color_pat.color_8, lw, lw, tipLength);
		break;
	case 9:
		cv::arrowedLine(img, p1, p2, color_pat.color_9, lw, lw, tipLength);
		break;
	case -1:
		cv::arrowedLine(img, p1, p2, color_pat.color_blue, lw, lw, tipLength);
		break;
	case -2:
		cv::arrowedLine(img, p1, p2, color_pat.color_green, lw, lw, tipLength);
		break;
	case -3:
		cv::arrowedLine(img, p1, p2, color_pat.color_red, lw, lw, tipLength);
		break;
	default:
		cv::arrowedLine(img, p1, p2, color_pat.color_default, lw, lw, tipLength);
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

Point_2 CCVWindow::convert_to_Point_2(cv::Point p)
{
	return Point_2(p.x, p.y);
}

