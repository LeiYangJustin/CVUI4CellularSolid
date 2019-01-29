#ifndef C_GEO_CALCULATOR_H
#define C_GEO_CALCULATOR_H

// Eigen
#include <Eigen/Core>
#include <Eigen/SVD>

// CGAL
#include "../DataColle/types.h"

//
#include "algprereq.h"


class ALGCOLLE_CLASS CGeoCalculator
{
public:
	CGeoCalculator();
	~CGeoCalculator();

	// cos angle
	static double compute_cos_angle(Point_2 p1, Point_2 p2);
	static double compute_cos_angle(Segment_2 s1, Segment_2 s2);
	static double compute_cos_angle(Vector_2 s1, Vector_2 s2);

	// sin angle
	static double compute_sin_angle(Vector_2 p1, Vector_2 p2);
	static double compute_sin_angle(Point_2 p1, Point_2 p2);

	// intersection of two lines
	static double intersection_of_two_lines(double x0, double y0, double x1, double y1,
		double x2, double y2, double x3, double y3);
	static double projection_of_point_to_line(Point_2 q, Point_2 p1, Point_2 p2);

	// interpolation
	template <class T>
	static double bilinear_interpl(Point_2 query, Point_2 lb_pt, Point_2 rt_pt,
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &F);

	template <class T>
	static double bilinear_interpl(Point_2 query, 
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &F,
		double maxX, double maxY, double minX = 0.0, double minY = 0.0);

	// transformation
	// rotation from S_right to S_left; R*S_right = S_left
	static Eigen::MatrixXd compute_rotation_from_right_to_left(Eigen::MatrixXd S0, Eigen::MatrixXd S1);
	static double compute_uniformScale_from_right_to_left(Eigen::MatrixXd S0, Eigen::MatrixXd S1);

	// eigen values
	static void compute_PCA_with_SVD(std::vector<Point_2> p_list,
		std::vector<std::vector<double>>& eigenvectors, std::vector<double>& eigenvalues);

	// CGAL
	// norm
	static double cgal_norm(Vector_2 v);
	static double cgal_norm(Point_2 v);



};

template<class T>
inline double CGeoCalculator::bilinear_interpl(Point_2 query, Point_2 lb_pt, Point_2 rt_pt, 
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& F)
{
	//https://en.wikipedia.org/wiki/Bilinear_interpolation
	// lb = left bottom, rt = right top
	double x1, x2, y1, y2, x, y;
	x1 = lb_pt.x(); x2 = rt_pt.x(); y1 = lb_pt.y(); y2 = rt_pt.y();
	x = query.x(); y = query.y();

	//
	Eigen::Matrix<T, 2, 1> vx;
	vx(0) = x2 - x; vx(1) = x - x1;
	Eigen::Matrix<T, 2, 1> vy; 
	vy(0) = y2 - y; vy(1) = y - y1;
	Eigen::Matrix<T, 2, 2> Mf;
	// four points Q11 = (x1, y1), Q12 = (x1, y2), Q21 = (x2, y1), and Q22 = (x2, y2).
	Mf << F(x1, y1), F(x1, y2),
		F(x2, y1), F(x2, y2);
	double val = vx.transpose()*Mf*vy;
	val /= (x2 - x1)*(y2 - y1);

	//
	return val;
}

template<class T>
inline double CGeoCalculator::bilinear_interpl(Point_2 query, 
	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& F, 
	double maxX, double maxY, double minX, double minY)
{
	double x1, x2, y1, y2, x, y;
	x1 = std::max(std::floor(query.x()), minX);
	y1 = std::max(std::floor(query.y()), minY);
	x2 = std::min(std::ceil(query.x()), maxX);
	y2 = std::min(std::ceil(query.y()), maxY);
	x = query.x(); 
	y = query.y();

	//std::cout << "x1, y1 " << x1 << ", " << y1 << std::endl;
	//std::cout << "x2, y2 " << x2 << ", " << y2 << std::endl;
	//std::cout << "x, y " << x << ", " << y << std::endl;

	//
	Eigen::Vector2d vx;
	vx(0) = x2 - x; vx(1) = x - x1;
	Eigen::Vector2d vy;
	vy(0) = y2 - y; vy(1) = y - y1;
	Eigen::Matrix2d Mf;

	double val = 0.0;
	if (vx.norm() == 0 || vy.norm() == 0)
		val = F(x1, y1);
	else if (vx(0) == 0 && vy(0) == 0)
		val = F(x1, y1);
	else if (vx(1) == 0 && vy(1) == 0)
		val = F(x2, y2);
	else {
		// four points Q11 = (x1, y1), Q12 = (x1, y2), Q21 = (x2, y1), and Q22 = (x2, y2).
		//std::cout << "bilinear interp" << std::endl;
		Mf << F(x1, y1), F(x1, y2),
			F(x2, y1), F(x2, y2);
		val = vx.transpose()*Mf*vy;
		val /= (x2 - x1)*(y2 - y1);
	}
	return val;
}

#endif