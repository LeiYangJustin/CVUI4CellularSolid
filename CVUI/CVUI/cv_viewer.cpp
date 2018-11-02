#include "cv_viewer.h"

void onMouseCViewer(int event, int x, int y, int flag, void* ptr)
{
	if (event != cv::EVENT_LBUTTONDOWN)
		return;

	std::cout << "click at" << x << " " << y << std::endl;
}

CViewer::CViewer()
{
}

CViewer::~CViewer()
{
}

void CViewer::draw(std::string winName, int elapseTime)
{
	cv::Mat tmp_img = data_.background_image.clone();
	//cv::imshow("aaa", tmp_img);
	//cv::waitKey(0);
	int lineWidth = 2;
	int colorID = 0;
	for (auto polyline : data_.polyline_list)
		drawPolyLine(tmp_img, polyline, data_.polyline_color_list[colorID++], lineWidth, true);

	colorID = 0;
	for (auto segment : data_.segment_list)
		drawSegment(tmp_img, segment.first, segment.second, data_.segment_color_list[colorID++], lineWidth);
	
	colorID = 0;
	for (auto pt : data_.single_pt_list)
		drawSegment(tmp_img, pt, pt, data_.single_pt_color_list[colorID++], lineWidth);

	colorID = 0;
	for (int i = 0; i < data_.ptset_list.size(); i++)
		drawPointSet(tmp_img, data_.ptset_list[i], data_.ptset_color_list[i], lineWidth);

	cv::namedWindow(winName);
	cv::imshow(winName, tmp_img);

	//cv::setMouseCallback(winName, onMouseCViewer);
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

	cv::waitKey(elapseTime);
	//cv::destroyWindow(winName);
}


void CViewer::drawSegment(cv::Mat & img, cv::Point p1, cv::Point p2, int label, int lw)
{
	COLOR_PALETTE color_pat;
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


void CViewer::drawPolyLine(cv::Mat & img, std::vector<cv::Point> plist, int label, int lw, bool is_closed)
{
	if (is_closed)
		plist.push_back(plist.front());

	COLOR_PALETTE color_pat;
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

void CViewer::drawPointSet(cv::Mat & img, std::vector<cv::Point> plist, std::vector<int> label_list, int lw)
{
	COLOR_PALETTE color_pat;
	float tipLength = 0.5;

	for (int i = 0; i < plist.size(); i++)
	{
		int label = label_list[i];
		switch (label % 10)
		{
		case 0:
			cv::line(img, plist[i], plist[i], color_pat.color_0, lw, lw, tipLength);
			break;
		case 1:
			cv::line(img, plist[i], plist[i], color_pat.color_1, lw, lw, tipLength);
			break;
		case 2:
			cv::line(img, plist[i], plist[i], color_pat.color_2, lw, lw, tipLength);
			break;
		case 3:
			cv::line(img, plist[i], plist[i], color_pat.color_3, lw, lw, tipLength);
			break;
		case 4:
			cv::line(img, plist[i], plist[i], color_pat.color_4, lw, lw, tipLength);
			break;
		case 5:
			cv::line(img, plist[i], plist[i], color_pat.color_5, lw, lw, tipLength);
			break;
		case 6:
			cv::line(img, plist[i], plist[i], color_pat.color_6, lw, lw, tipLength);
			break;
		case 7:
			cv::line(img, plist[i], plist[i], color_pat.color_7, lw, lw, tipLength);
			break;
		case 8:
			cv::line(img, plist[i], plist[i], color_pat.color_8, lw, lw, tipLength);
			break;
		case 9:
			cv::line(img, plist[i], plist[i], color_pat.color_9, lw, lw, tipLength);
			break;
		case -1:
			cv::line(img, plist[i], plist[i], color_pat.color_blue, lw, lw, tipLength);
			break;
		case -2:
			cv::line(img, plist[i], plist[i], color_pat.color_green, lw, lw, tipLength);
			break;
		case -3:
			cv::line(img, plist[i], plist[i], color_pat.color_red, lw, lw, tipLength);
			break;
		default:
			cv::line(img, plist[i], plist[i], color_pat.color_default, lw, lw, tipLength);
			break;
		}
	}
	
}
