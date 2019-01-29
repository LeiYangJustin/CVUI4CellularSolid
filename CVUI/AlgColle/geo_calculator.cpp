#include "geo_calculator.h"



CGeoCalculator::CGeoCalculator()
{
}


CGeoCalculator::~CGeoCalculator()
{
}

double CGeoCalculator::compute_cos_angle(Point_2 p1, Point_2 p2)
{	
	return compute_cos_angle(p1 - Point_2(0, 0), p2 - Point_2(0, 0));
}

double CGeoCalculator::compute_cos_angle(Segment_2 s1, Segment_2 s2)
{
	return compute_cos_angle(s1.to_vector(), s2.to_vector());
}

double CGeoCalculator::compute_cos_angle(Vector_2 s1, Vector_2 s2)
{
	return s1 *s2 / sqrt(s1.squared_length()*s2.squared_length());
}

double CGeoCalculator::compute_sin_angle(Vector_2 v1, Vector_2 v2)
{
	v1 /= sqrt(v1.squared_length());
	v2 /= sqrt(v1.squared_length());
	return (v1.x()*v2.y() - v1.y()*v2.x()) / 2.0;
}

double CGeoCalculator::compute_sin_angle(Point_2 p1, Point_2 p2)
{
	Vector_2 v1 = p1 - Point_2(0, 0);
	Vector_2 v2 = p1 - Point_2(0, 0);
	v1 /= sqrt(v1.squared_length());
	v2 /= sqrt(v1.squared_length());
	return (v1.x()*v2.y() - v1.y()*v2.x()) / 2.0;
}

double CGeoCalculator::intersection_of_two_lines(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3)
{
	Eigen::Matrix2d A;
	A << x1 - x0, x2 - x3,
		y1 - y0, y2 - y3;
	Eigen::Vector2d b;
	b << x2 - x0,
		y2 - y0;
	Eigen::Vector2d t = A.colPivHouseholderQr().solve(b);
	double tp = t(1);
	return tp;
}

double CGeoCalculator::projection_of_point_to_line(Point_2 q, Point_2 p1, Point_2 p2)
{
	Vector_2 v1 = q - p1;
	Vector_2 v2 = p2 - p1;
	return v1*v2 / sqrt(v2.squared_length());
}

Eigen::MatrixXd CGeoCalculator::compute_rotation_from_right_to_left(Eigen::MatrixXd S0, Eigen::MatrixXd S1)
{
	assert(S1.cols() == S0.cols());
	// assuming S0 and S1 are already shifted to their centers (means)
	// we want to compute a rotation R that rotates S1 (A^T) to S0 (B^T)
	// | RA - B |_{Frobenius_norm} -> min, where R is an orthogonal matrix
	//https://www.ltu.se/cms_fs/1.51590!/svd-fitting.pdf

	Eigen::MatrixXd C = S0.transpose()*S1;
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(C, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::MatrixXd R = svd.matrixU()*svd.matrixV().transpose();
	return R;
}

double CGeoCalculator::compute_uniformScale_from_right_to_left(Eigen::MatrixXd S0, Eigen::MatrixXd S1)
{
	assert(S0.cols() == S1.cols());
	
	// svd
	Eigen::JacobiSVD<Eigen::MatrixXd> svd_s0(S0, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::JacobiSVD<Eigen::MatrixXd> svd_s1(S1, Eigen::ComputeFullU | Eigen::ComputeFullV);
	
	// compute averaged scale
	double s = 0.0;
	int dim = S0.cols();
	for (int i = 0; i < dim; i++) {
		// Eigen::JacobiSVD<>::singularValues is a vector
		s += svd_s0.singularValues()(i)/svd_s1.singularValues()(i);
	}
	s /= dim;
	return s;
}


double CGeoCalculator::cgal_norm(Vector_2 v)
{
	return sqrt(v.x()*v.x() + v.y()*v.y());
}

double CGeoCalculator::cgal_norm(Point_2 v)
{
	return sqrt(v.x()*v.x() + v.y()*v.y());
}

void CGeoCalculator::compute_PCA_with_SVD(std::vector<Point_2> p_list, std::vector<std::vector<double>>& eigenvectors, std::vector<double>& eigenvalues)
{
	eigenvectors.clear();
	eigenvalues.clear();

	int Dim = 2;
	Eigen::MatrixXd P(p_list.size(), Dim);
	Eigen::MatrixXd V(Dim, Dim);
	Eigen::VectorXd S(Dim);

	Vector_2 center_p(0, 0, 0);
	for (int i = 0; i < p_list.size(); i++)
	{
		center_p += (p_list[i]-Point_2(0, 0));
	}
	center_p = center_p / double(p_list.size());
	for (int i = 0; i < p_list.size(); i++)
	{
		P(i, 0) = p_list[i].x() - center_p.x();
		P(i, 1) = p_list[i].y() - center_p.y();
		//P(i, 2) = p_list[i][2] - center_p[2];
	}
	Eigen::MatrixXd CovP = P.transpose() * P;
	CovP = CovP / double(p_list.size());
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(CovP, Eigen::ComputeThinU | Eigen::ComputeThinV);
	V = svd.matrixV();
	S = svd.singularValues();

	//std::cout << "using ComputePCAwithSVD function" << std::endl;
	//std::cout << "S = [" << S << "]" << std::endl;
	//std::cout << "V = [" << V << "]" << std::endl;

	for (int i = 0; i < Dim; i++)
	{
		std::vector<double> tmp_vector;
		auto n = V.col(i);
		n.normalized();
		tmp_vector.push_back(n[0]);
		tmp_vector.push_back(n[1]);
		//tmp_vector.push_back(n[2]);
		eigenvectors.push_back(tmp_vector);
		double s = sqrt(S(i));
		eigenvalues.push_back(s);
	}
}
