#include "../DataColle/voronoi_diagram.h"
#include <vector>

// INCLUDE OPENCV DIR
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include <opencv2/opencv.hpp>

class CVoronoiDrawer
{
public:
	CVoronoiDrawer();
	~CVoronoiDrawer();
	std::vector<std::pair<cv::Point, cv::Point>> DrawVoronoi();
	void SetVD(CVoronoiDiagram* vd);

private:
	CVoronoiDiagram *vd_;
	cv::Point convert_to_cvPoint(Point2);
};