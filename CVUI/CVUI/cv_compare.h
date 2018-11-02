#ifndef CV_COMPARE_H
#define CV_COMPARE_H

// INCLUDE OPENCV DIR
#include <opencv2/core/core.hpp>

struct CV_LESS_COMPARE
{
	bool operator() (const cv::Point & p1, const cv::Point & p2) const
	{
		if (p1.x < p2.x)
			return true;
		else if (p1.x > p2.x)
			return false;
		else
		{
			if (!(p1.y > p2.y))
				return true;
			else
				return false;
		}
	}
};

#endif // !CV_LESS_COMPARE_H
