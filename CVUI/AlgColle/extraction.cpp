#include "extraction.h"

#include <iostream>
#include <fstream>
#include <iterator>     // std::istream_iterator

void CExtraction::setup_background(int rows, int cols, std::vector<double> sz, std::vector<double> vz)
{
	std::cout << "rows = " << rows << ", cols = " << cols << std::endl;

	rows_ = rows;
	cols_ = cols;

	x_grid_.clear();
	y_grid_.clear();

	for (int i = 0; i < cols_; i++)
	{
		x_grid_.push_back(i);
	}
	for (int i = 0; i < rows_; i++)
	{
		y_grid_.push_back(i);
	}

	assert(sz.size() == cols_*rows_);
	assert(sz.size() == vz.size());

	Eigen::MatrixXd solid_xy_z = Eigen::Map<Eigen::MatrixXd>(sz.data(), rows_, cols_);
	Eigen::MatrixXd void_xy_z = Eigen::Map<Eigen::MatrixXd>(vz.data(), rows_, cols_);
	z_ = void_xy_z - solid_xy_z;

	compute_field_gradients(z_, zgx_, zgy_);

	//CMeshIO::write_input_fields_for_check(z_, zgx_, zgy_);
}

void CExtraction::run_optimize()
{
	// MAX ITER FOR THE GLOBAL LOOP
	const int MAX_ITERS = 100;

	// PARAMETERS FOR SYNTHESIS
	int verbose = 0;
	int stepX = 0.0;
	int stepW = 0.0;
	int epsilon = 1.0;
	int frequency = 0;
	int max_newton_iters = 500;
	int max_opt_iters = 500;

	// SET THE DOMAIN
	Domain in_domain;
	in_domain.init_rectangle(cols_, rows_);
	in_domain.init_area();
	scene_.set_domain(in_domain);

	// INIT THE SITES
	int nb = 100;
	scene_.generate_random_sites(nb);

	//
	int ext_iter = 0;
	while (ext_iter++ < MAX_ITERS)
	{
		// SCENE_OPTIMIZE
		scene_.optimize_all_with_background_constraint(stepW, stepX, 
			max_newton_iters, 
			epsilon, 
			max_opt_iters, 
			std::cout);

		scene_.get_rt(rt_);

		// INSERT PTS
		remove_seeds_in_proximity();
		insert_seeds_at_edges();
		insert_seeds_at_faces();
		scene_.set_rt(rt_);
	}

	std::cout << "finished extraction..." << std::endl;

}

void CExtraction::remove_seeds_in_proximity(double merging_length)
{
	// REMOVE
	std::vector<RT::Weighted_point> compact_seeds;
	for (auto vit = rt_->finite_vertices_begin(); vit != rt_->finite_vertices_end(); ++vit)
	{
		bool has_vertex_in_proximity = false;
		double hpt = dist_transform_field_S(Eigen::Vector2d(vit->point().x(), vit->point().y()));
		//WPoint wpt = WPoint(pt, 0);
		for (auto seed : compact_seeds)
		{
			double dist = sqrt((vit->point().point() - seed.point()).squared_length());
			if (dist < merging_length)
			{
				double hs = dist_transform_field_S(Eigen::Vector2d(seed.point().x(), seed.point().y()));
				if (hs < hpt)
				{
					seed = vit->point();
					has_vertex_in_proximity = true;
					break;
				}
			}
		}
		if (!has_vertex_in_proximity) {
			compact_seeds.push_back(vit->point());
		}
	}

	rt_->clear();
	rt_->insert(compact_seeds.begin(), compact_seeds.end());
}

void CExtraction::insert_seeds_at_edges()
{
	std::vector<Weighted_point_2> seeds_to_insert;

	// for each edge in the new built triangulation, 
	// if the edge passes a local minimum, then we shall add a point there to describe it

	for (auto eit = rt_->all_edges_begin(); eit != rt_->all_edges_end(); ++eit)
	{
		// skip infinite edges
		if (rt_->is_infinite(eit))
			continue;

		auto fh = eit->first;
		auto vh = fh->vertex(eit->second);
		auto psrc = fh->vertex(fh->ccw(eit->second))->point().point();
		auto ptgt = fh->vertex(fh->cw(eit->second))->point().point();

		// sampling on edge(psrc, ptgt)
		Vector_2 vdir = ptgt - psrc;
		if (vdir.squared_length() < 0.00001)
		{
			// should not have these edges because of merging is done ahead of insersion
			// skip very short edges
			std::cerr << "very short edge" << std::endl;
			//continue;
		}
		// we only sample those in the middle, because the two end should have grad = 0
		double step = 0.1;
		double min_grad_q_norm = 4;
		double min_t = -1;
		int cnt_turning_sign = 0;
		double sign = -1;
		double h0 = dist_transform_field_S(Eigen::Vector2d(psrc.x(), psrc.y()));
		for (double t = 0.05; t < 0.95; t += step)
		{
			Point_2 q = psrc + t*vdir;
			double hq = dist_transform_field_S(Eigen::Vector2d(q.x(), q.y()));
			if (sign * hq < 0)
			{
				sign = -sign;
				cnt_turning_sign++;
			}
			Eigen::Vector2d grad_q = dist_transform_grad_S(Eigen::Vector2d(q.x(), q.y()));
			if (hq < 0 && grad_q.norm() < min_grad_q_norm && cnt_turning_sign == 2)
			{
				min_t = t;
				min_grad_q_norm = grad_q.norm();
			}
		}
		if (min_t > 0 && cnt_turning_sign > 3)
		{
			auto p = psrc + min_t*vdir;
			double hinsert = dist_transform_field_S(Eigen::Vector2d(p.x(), p.y()));
			seeds_to_insert.push_back(Weighted_point_2(p, 0));
		}
	}

	std::cout << "# new inserted seeds: " << seeds_to_insert.size() << std::endl;

	// INSERT
	rt_->insert(seeds_to_insert.begin(), seeds_to_insert.end());

}

void CExtraction::insert_seeds_at_faces()
{
	std::cout << "under construction" << std::endl;
}

void CExtraction::compute_field_gradients(
	const Eigen::MatrixXd &zmap,
	Eigen::MatrixXd &zgradx,
	Eigen::MatrixXd &zgrady)
{
	zgradx.resize(zmap.rows(), zmap.cols());
	zgradx.setZero();

	zgrady.resize(zmap.rows(), zmap.cols());
	zgrady.setZero();

	double gx, gy;
	for (int ix = 1; ix + 1 < x_grid_.size(); ix++)
	{
		for (int iy = 1; iy + 1 < y_grid_.size(); iy++)
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

double CExtraction::get_interp_val_from_2dgrid(const Eigen::Vector2d & pos, const Eigen::MatrixXd & map)
{
	double x, y;
	x = pos[0];
	y = pos[1];

	double lx, ux, ly, uy;
	int idx, idy;
	idx = -1; idy = -1;

	// X
	for (int i = 0; i + 1 < x_grid_.size(); i++)
	{
		//
		if (x < x_grid_[0])
		{
			idx = 0;
			lx = x_grid_[0];
			ux = x_grid_[1];
		}
		//
		if (x_grid_[i] <= x && x < x_grid_[i + 1])
		{
			idx = i;
			lx = x_grid_[i];
			ux = x_grid_[i + 1];
		}
		// 
		if (x >= x_grid_.back())
		{
			idx = x_grid_.size() - 2;
			lx = x_grid_.size() - 2;
			ux = x_grid_.size() - 1;
		}
	}


	// Y
	for (int i = 0; i + 1 < y_grid_.size(); i++)
	{
		//
		if (y < y_grid_[0])
		{
			idy = 0;
			ly = y_grid_[0];
			uy = y_grid_[1];
		}
		//
		if (y_grid_[i] <= y && y < y_grid_[i + 1])
		{
			idy = i;
			ly = y_grid_[i];
			uy = y_grid_[i + 1];
		}
		// 
		if (y >= y_grid_.back())
		{
			idy = y_grid_.size() - 2;
			ly = y_grid_.size() - 2;
			uy = y_grid_.size() - 1;
		}
	}

	assert(idx > -1);
	assert(idy > -1);
	// bilinear interpolation
	double s = (x - lx) / (ux - lx);
	double t = (y - ly) / (uy - ly);

	// check this: x -> cols, 
	int stride = x_grid_.size();
	double v_0 = map(idy, idx);
	double v_1 = map(idy + 1, idx);
	double v_2 = map(idy + 1, idx + 1);
	double v_3 = map(idy, idx + 1);

	double v_interp = (1 - s)*(1 - t)*v_0 + s*(1 - t)*v_1 + s*t*v_2 + (1 - s)*t*v_3;

	return v_interp;
}

