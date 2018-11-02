#ifndef C_VIEWER_H
#define C_VIEWER_H

#include "../DataColle/cgal_def.h"
#include "../DataColle/img_data.h"
#include <vector>
#include <iostream>

// INCLUDE OPENCV DIR
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include <opencv2/opencv.hpp>
//
//struct CV_LESS_COMPARE
//{
//	bool operator() (const cv::Point & p1, const cv::Point & p2) const
//	{
//		if (p1.x < p2.x)
//			return true;
//		else if (p1.x > p2.x)
//			return false;
//		else
//		{
//			if (!(p1.y > p2.y))
//				return true;
//			else
//				return false;
//		}
//	}
//};

struct COLOR_PALETTE {
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
	cv::Scalar color_blue = cv::Scalar(255, 0, 0);
	cv::Scalar color_green = cv::Scalar(0, 255, 0);
	cv::Scalar color_red = cv::Scalar(0, 0, 255);
};

struct ViewerContent
{
	ViewerContent() :
		background_image(cv::Mat(cv::Size(500, 500), CV_8UC3, cv::Scalar(255, 255, 255)))
	{
	};
	ViewerContent(int w, int h) :
		background_image(cv::Mat(cv::Size(w, h), CV_8UC3, cv::Scalar(255, 255, 255)))
	{
	};

	void clearSegmentList() {
		segment_color_list.clear();
		segment_list.clear();
	}
	void clearPolylineList() {
		polyline_color_list.clear();
		polyline_list.clear();
	}
	void clearSinglePtList() {
		single_pt_color_list.clear();
		single_pt_list.clear();
	}
	void clearPtSetList() {
		ptset_color_list.clear();
		ptset_list.clear();
	}

	cv::Mat background_image;
	std::vector< std::vector<cv::Point> > polyline_list;
	std::vector< int > polyline_color_list;
	std::vector< std::pair<cv::Point, cv::Point> > segment_list;
	std::vector< int > segment_color_list;
	std::vector< cv::Point > single_pt_list;
	std::vector< int > single_pt_color_list;
	std::vector< std::vector<cv::Point> > ptset_list;
	std::vector< std::vector<int> > ptset_color_list;
};

class CViewer
{
public:
	CViewer();
	CViewer(CImgData img_data){ data_.background_image = img_data.GetBackgroundImage(); }
	CViewer(int w, int h) :data_(w, h) {};
	~CViewer();
	
	void setImgData(CImgData img_data)
	{
		data_.background_image = img_data.GetSrcImg();
	}

	int addPolyLine(std::vector<Point_2> polyline, int colorID = -1)
	{
		std::vector<cv::Point> cv_polyline;
		for (auto pt : polyline)
			cv_polyline.push_back(cv::Point(pt.x(), pt.y()));
		data_.polyline_list.push_back(cv_polyline);
		data_.polyline_color_list.push_back(colorID);
		return (data_.polyline_list.size() - 1);
	}

	int addSegment(std::pair<Point_2, Point_2> segment, int colorID = -2)
	{
		cv::Point cvp1(segment.first.x(), segment.first.y());
		cv::Point cvp2(segment.second.x(), segment.second.y());
		data_.segment_list.push_back(std::make_pair(cvp1, cvp2));
		data_.segment_color_list.push_back(colorID);
		return (data_.segment_list.size() - 1);
	}
	int addSegment(Segment_2 segment, int colorID = -2)
	{
		cv::Point cvp1(segment.source().x(), segment.source().y());
		cv::Point cvp2(segment.target().x(), segment.target().y());
		data_.segment_list.push_back(std::make_pair(cvp1, cvp2));
		data_.segment_color_list.push_back(colorID);
		return (data_.segment_list.size() - 1);
	}

	int addPoint(Point_2 pt, int colorID = -3)
	{
		data_.single_pt_list.push_back(cv::Point(pt.x(), pt.y()));
		data_.single_pt_color_list.push_back(colorID);
		return(data_.single_pt_list.size() - 1);
	}

	int addPointSet(std::vector<Point_2> pts, std::vector<int> colorIDs)
	{
		std::vector<cv::Point> cvpts;
		for(auto pt : pts)
			cvpts.push_back(cv::Point(pt.x(), pt.y()));
		data_.ptset_list.push_back(cvpts);
		data_.ptset_color_list.push_back(colorIDs);
		return(data_.ptset_list.size() - 1);
	}
	
	void draw(std::string winName, int elapseTime = 0);

	void clearData()
	{
		data_.clearPolylineList();
		data_.clearPtSetList();
		data_.clearSegmentList();
		data_.clearSinglePtList();
	}

private:
	ViewerContent data_;

	void drawSegment(cv::Mat & img, cv::Point p1, cv::Point p2, int label, int lw);
	void drawPolyLine(cv::Mat & img, std::vector<cv::Point> plist, int label, int lw, bool is_closed);
	void drawPointSet(cv::Mat & img, std::vector<cv::Point> plist, std::vector<int> label_list, int lw);


};

#endif // !C_VIEWER_H


