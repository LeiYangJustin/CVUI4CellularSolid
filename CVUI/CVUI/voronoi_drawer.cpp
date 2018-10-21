#include "voronoi_drawer.h"

CVoronoiDrawer::CVoronoiDrawer()
{
	vd_ = new CVoronoiDiagram;
}

CVoronoiDrawer::~CVoronoiDrawer()
{
	delete vd_;
}

std::vector<std::pair<cv::Point, cv::Point>> CVoronoiDrawer::DrawVoronoi()
{
	std::vector<std::pair<Point2, Point2>> voronoi_segments;
	vd_->GetCroppedVoronoiSegments(voronoi_segments);

	std::vector<std::pair<cv::Point, cv::Point>> cv_voronoi_segments;
	for (int i = 0; i < voronoi_segments.size(); i++)
	{
		cv::Point cvp1 = convert_to_cvPoint(voronoi_segments[i].first);
		cv::Point cvp2 = convert_to_cvPoint(voronoi_segments[i].second);
		cv_voronoi_segments.push_back(std::make_pair(cvp1, cvp2));
	}
	return cv_voronoi_segments;
}

void CVoronoiDrawer::SetVD(CVoronoiDiagram * vd)
{
	vd_ = vd;
}

cv::Point CVoronoiDrawer::convert_to_cvPoint(Point2 p)
{
	return cv::Point(p.x(), p.y());
}
