#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"
#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>
using namespace std;
using namespace cv;
using namespace std::chrono;
// Function declarations
void drawAxis(Mat&, Point, Point, Scalar, const float);
double getOrientation(const vector<Point> &, Mat&);
void drawAxis(Mat& img, Point p, Point q, Scalar colour, const float scale = 0.2)
{
	double angle = atan2((double)p.y - q.y, (double)p.x - q.x); // angle in radians
	double hypotenuse = sqrt((double)(p.y - q.y) * (p.y - q.y) + (p.x - q.x) * (p.x - q.x));
	// Here we lengthen the arrow by a factor of scale
	q.x = (int)(p.x - scale * hypotenuse * cos(angle));
	q.y = (int)(p.y - scale * hypotenuse * sin(angle));
	line(img, p, q, colour, 1, LINE_AA);
	// create the arrow hooks
	p.x = (int)(q.x + 9 * cos(angle + CV_PI / 4));
	p.y = (int)(q.y + 9 * sin(angle + CV_PI / 4));
	line(img, p, q, colour, 1, LINE_AA);
	p.x = (int)(q.x + 9 * cos(angle - CV_PI / 4));
	p.y = (int)(q.y + 9 * sin(angle - CV_PI / 4));
	line(img, p, q, colour, 1, LINE_AA);
}
double getOrientation(const vector<Point> &pts, Mat &img)
{
	//Construct a buffer used by the pca analysis
	int sz = static_cast<int>(pts.size());
	Mat data_pts = Mat(sz, 2, CV_64F);
	for (int i = 0; i < data_pts.rows; i++)
	{
		data_pts.at<double>(i, 0) = pts[i].x;
		data_pts.at<double>(i, 1) = pts[i].y;
	}
	//Perform PCA analysis
	PCA pca_analysis(data_pts, Mat(), PCA::DATA_AS_ROW);
	//Store the center of the object
	Point cntr = Point(static_cast<int>(pca_analysis.mean.at<double>(0, 0)),
		static_cast<int>(pca_analysis.mean.at<double>(0, 1)));
	//Store the eigenvalues and eigenvectors
	vector<Point2d> eigen_vecs(2);
	vector<double> eigen_val(2);
	for (int i = 0; i < 2; i++)
	{
		eigen_vecs[i] = Point2d(pca_analysis.eigenvectors.at<double>(i, 0),
			pca_analysis.eigenvectors.at<double>(i, 1));
		eigen_val[i] = pca_analysis.eigenvalues.at<double>(i);
		
		//std::cout << "eigen_val[" << i << "] = " << eigen_val[i] << std::endl;
	}

	// Draw the principal components
	circle(img, cntr, 3, Scalar(255, 0, 255), 2);
	Point p1 = cntr + 0.2 * Point(static_cast<int>(eigen_vecs[0].x * eigen_val[0]), static_cast<int>(eigen_vecs[0].y * eigen_val[0]));
	Point p2 = cntr - 0.2 * Point(static_cast<int>(eigen_vecs[1].x * eigen_val[1]), static_cast<int>(eigen_vecs[1].y * eigen_val[1]));
	drawAxis(img, cntr, p1, Scalar(0, 255, 0), 1);
	drawAxis(img, cntr, p2, Scalar(255, 255, 0), 1);
	double angle = atan2(eigen_vecs[0].y, eigen_vecs[0].x); // orientation in radians
	return angle;
}
int main(int argc, char** argv)
{
	// Load image
	std::string filename = "D:\\MyProjects\\CVUI4CellularSolid\\CVUI\\img_data\\example2.png";
	Mat src = imread(filename);
	// Check if image is loaded successfully
	if (src.empty())
	{
		cout << "Problem loading image!!!" << endl;
		return EXIT_FAILURE;
	}
	imshow("src", src);
	// Convert image to grayscale
	Mat gray;
	cvtColor(src, gray, COLOR_BGR2GRAY);
	// Convert image to binary
	Mat bw;
	threshold(gray, bw, 50, 255, THRESH_BINARY | THRESH_OTSU);

	// pick a point
	cv::Point upperleft_crn(250, 250);
	int width = 101;
	int height = 101;
	cv::Point cntr = upperleft_crn + cv::Point(width / 2, height / 2);

	//high_resolution_clock::time_point t1 = high_resolution_clock::now();

	for (int k = 0; k < 50000; ++k)
	{
		// crop the image centered at pt with width and height as defined
		Mat bw_patch;
		getRectSubPix(bw, cv::Size(width, height), cntr, bw_patch);
		bitwise_not(bw_patch, bw_patch);
		// get contours
		vector<vector<Point> > contours;
		findContours(bw_patch, contours, RETR_EXTERNAL, CHAIN_APPROX_NONE);
		for (size_t i = 0; i < contours.size(); i++)
		{
			// shift contours from bw_patch to src
			for (size_t j = 0; j < contours[i].size(); j++)
			{
				contours[i][j] += upperleft_crn;
			}
			double in_contour = pointPolygonTest(contours[i], cntr, true);
			//std::cout << in_contour << std::endl;
			if (in_contour <= 0)
				continue;
			//// Draw each contour only for visualisation purposes
			//drawContours(src, contours, static_cast<int>(i), Scalar(0, 0, 255), 2);
			// Find the orientation of each shape
			getOrientation(contours[i], src);

			//high_resolution_clock::time_point t1 = high_resolution_clock::now();
			//high_resolution_clock::time_point t2 = high_resolution_clock::now();
			//duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
			//std::cout << "It took me " << time_span.count() << " seconds.";
			//std::cout << std::endl;
		}
	}
	//high_resolution_clock::time_point t2 = high_resolution_clock::now();
	//duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	//std::cout << "It took me " << time_span.count() << " seconds.";
	//std::cout << std::endl;
	circle(src, cntr, 5, Scalar(0, 255, 0), 2);
	imshow("output", src);
	waitKey();
	return 0;
}