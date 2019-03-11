#include "extraction.h"

#include <iostream>
#include <fstream>
#include <iterator>     // std::istream_iterator
#include <ctime>
#include <ratio>
#include <chrono>

//void CMeshExtractor::setup_background(int rows, int cols, std::vector<double> &sz, std::vector<double> &vz)
//{
//	std::cout << "rows = " << rows << ", cols = " << cols << std::endl;
//
//	rows_ = rows;
//	cols_ = cols;
//
//	assert(sz.size() == cols_*rows_);
//	assert(sz.size() == vz.size());
//
//	Eigen::MatrixXd solid_xy_z = Eigen::Map<Eigen::MatrixXd>(sz.data(), rows_, cols_);
//	Eigen::MatrixXd void_xy_z = Eigen::Map<Eigen::MatrixXd>(vz.data(), rows_, cols_);
//	z_ = void_xy_z - solid_xy_z;
//
//	compute_field_gradients(z_, zgy_, zgx_);
//
//	bg_scene_.construct_bg_info(cols_, rows_, z_, zgx_, zgy_);
//	
//	write_input_fields_for_check(z_, zgx_, zgy_);
//}

void CMeshExtractor::setup_background(int rows, int cols, 
	std::vector<double>& sz, std::vector<double>& vz, std::vector<std::vector<double>> in_edgePts)
{
	std::cout << "rows = " << rows << ", cols = " << cols << std::endl;

	rows_ = rows;
	cols_ = cols;

	assert(sz.size() == cols_*rows_);
	assert(sz.size() == vz.size());

	// distance field
	Eigen::MatrixXd solid_xy_z = Eigen::Map<Eigen::MatrixXd>(sz.data(), rows_, cols_);
	Eigen::MatrixXd void_xy_z = Eigen::Map<Eigen::MatrixXd>(vz.data(), rows_, cols_);
	z_ = void_xy_z - solid_xy_z;

	// gradient field
	compute_field_gradients(z_, zgy_, zgx_);
	
	// set info to bg_scene
	std::vector<Point_2> edgePts;
	for (int i = 0; i < in_edgePts.size(); ++i)
	{
		// convert a cv::Point format to a Point_2, switching row-col to x-y
		edgePts.push_back(Point_2(in_edgePts[i][1], in_edgePts[i][0]));
	}
	bg_scene_.construct_bg_info(cols_, rows_, z_, zgx_, zgy_, edgePts);
	write_input_fields_for_check(z_, zgx_, zgy_);
}

void CMeshExtractor::run_extraction()
{
	// PARAMETERS FOR SYNTHESIS
	int verbose = 0;
	double stepX = 0.0;
	double stepW = 0.0;
	double epsilon = 1.0;
	/*int frequency = 0;*/
	int max_newton_iters = 5;
	int max_opt_iters = 100;

	// SET THE DOMAIN
	bg_scene_.set_domain();

	// INIT THE SITES
	int nb = 30;
	bg_scene_.generate_random_sites(nb);

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	// SCENE_OPTIMIZE
	bg_scene_.bg_optimize_all_for_extraction(stepW, stepX,
		max_newton_iters, 
		epsilon, 
		max_opt_iters, 
		std::cout);

	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	std::cout << "It took me " << time_span.count() << " seconds.";
	std::cout << std::endl;

	std::cout << "finished extraction..." << std::endl;

}

bool CMeshExtractor::write_input_fields_for_check(Eigen::MatrixXd & z, Eigen::MatrixXd & zgx, Eigen::MatrixXd & zgy)
{
	std::cout << "write solid field" << std::endl;
	std::ofstream out_sz_file;
	out_sz_file.open("results\\solid_field.array");
	if (out_sz_file.is_open())
	{
		out_sz_file << z << std::endl;
		out_sz_file.close();
	}
	else {
		std::cout << "wrong at writing solid field" << std::endl;
		return false;
	}

	std::cout << "write solid grad X field" << std::endl;
	std::ofstream out_sgx_file;
	out_sgx_file.open("results\\solid_grad_x.array");
	if (out_sgx_file.is_open())
	{
		out_sgx_file << zgx << std::endl;
		out_sgx_file.close();
	}
	else {
		std::cout << "wrong at writing solid grad X field" << std::endl;
		return false;
	}

	std::cout << "write solid grad Y field" << std::endl;
	std::ofstream out_sgy_file;
	out_sgy_file.open("results\\solid_grad_y.array");
	if (out_sgy_file.is_open())
	{
		out_sgy_file << zgy << std::endl;
		out_sgy_file.close();
	}
	else {
		std::cout << "wrong at writing solid grad Y field" << std::endl;
		return false;
	}
	return true;
}

void CMeshExtractor::compute_field_gradients(
	const Eigen::MatrixXd &zmap,
	Eigen::MatrixXd &zgradx,
	Eigen::MatrixXd &zgrady)
{
	zgradx.resize(zmap.rows(), zmap.cols());
	zgradx.setZero();

	zgrady.resize(zmap.rows(), zmap.cols());
	zgrady.setZero();

	std::vector<double> x_grid;
	std::vector<double> y_grid;
	for (int i = 0; i < cols_; i++)
	{
		x_grid.push_back(i);
	}
	for (int i = 0; i < rows_; i++)
	{
		y_grid.push_back(i);
	}

	double gx, gy;
	for (int ix = 1; ix + 1 < x_grid.size(); ix++)
	{
		for (int iy = 1; iy + 1 < y_grid.size(); iy++)
		{
			// central difference
			gy = (zmap(iy + 1, ix - 1) - zmap(iy - 1, ix - 1)
				+ zmap(iy + 1, ix) - zmap(iy - 1, ix)
				+ zmap(iy + 1, ix + 1) - zmap(iy - 1, ix + 1)) / 3.0;

			gx = (zmap(iy - 1, ix + 1) - zmap(iy - 1, ix - 1)
				+ zmap(iy, ix + 1) - zmap(iy, ix - 1)
				+ zmap(iy + 1, ix + 1) - zmap(iy + 1, ix - 1)) / 3.0;

			zgradx(iy, ix) = gx;
			zgrady(iy, ix) = gy;
		}
	}
}

void CMeshExtractor::compute_principal_directions(const Eigen::MatrixXd & zmap, 
	Eigen::MatrixXd & pc1x, Eigen::MatrixXd & pc1y, 
	Eigen::MatrixXd & pc2x, Eigen::MatrixXd & pc2y)
{
	
	
}


