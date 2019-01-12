#include "mesh_optimizer.h"

#include <iostream>
#include <fstream>
#include <iterator>     // std::istream_iterator

void CMeshOptimizer::background_setup(int rows, int cols, std::vector<double> sz, std::vector<double> vz)
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

	//z_.resize(rows_, cols_);
	//z_.setZero();

	//for (int i = 0; i < z_.rows(); i++) {
	//	for (int j = 0; j < z_.cols(); j++)
	//	{
	//		z_(i, j) = z_(i, j) * z_(i, j) * z_(i, j);
	//	}
	//}


	compute_field_gradients(z_, zgx_, zgy_);

	write_input_fields_for_check();
}

void CMeshOptimizer::run_optimize()
{
	/* PARAMETERS */
	// weight for valley constraint
	const double alpha = 0.5;
	// weight for ridge constraint
	const double beta = 0.01;
	// weight for ODT part
	const double gamma_optimize_odt = 10;
	const double gamma_optimize_cvt = 1; // do not change this parameter for cvt (default = 1)
	const double gamma_refine = 0.01;
	// initial gradient step
	const double init_stepsize = 5;
	// length threshold to merge two sites
	const double merging_length = 10;
	// grid sampling spacing
	const double sampling_length = 10;
	// max iteration
	const int max_iter = 20;
	// refinement  parameters
	const double refine_thres = 0.005; // set to negative if no refinement is needed;
	const int max_refine_iter = 2;
	// outer bounding box
	const double offset = 0;
	// use CVT or ODT energy
	const bool use_cvt = false;

	/* SAMPLING THE DOMAIN */
	// bounding box
	std::vector<std::pair<Point_2, Point_2>> bounding_box;
	bounding_box.push_back(std::make_pair(Point_2(0,         0        ), Point_2(0,         rows_ - 1)));
	bounding_box.push_back(std::make_pair(Point_2(0,         rows_ - 1), Point_2(cols_ - 1, rows_ - 1)));
	bounding_box.push_back(std::make_pair(Point_2(cols_ - 1, rows_ - 1), Point_2(cols_ - 1, 0        )));
	bounding_box.push_back(std::make_pair(Point_2(cols_ - 1, 0),         Point_2(0,         0        )));
	// boundary sampling
	std::vector<Point_2> seeds;
	for (auto bound : bounding_box)
	{
		std::vector<Point_2> seeds_on_bound;
		boundary_sampling(bound.first, bound.second, seeds_on_bound);
		seeds.insert(seeds.end(), seeds_on_bound.begin(), seeds_on_bound.end());
	}
	// grid sampling
	for (int ix = 20; ix < cols_ - 20; ix += sampling_length)
	{
		for (int iy = 20; iy < rows_-20; iy += sampling_length)
		{
			seeds.push_back(Point_2(ix, iy));
		}
	}
	// outer bounding box; should be fixed
	seeds.push_back(Point_2(0 - offset, 0 - offset));
	seeds.push_back(Point_2(0 - offset, rows_ - 1 + offset));
	seeds.push_back(Point_2(cols_ - 1 + offset, rows_ - 1 + offset));
	seeds.push_back(Point_2(cols_ - 1 + offset, 0 - offset));

	/* MAIN SECTION FOR CONSTRAINED ODT OPTIMIZATION */
	double gamma_optimize = 0;
	if (use_cvt)
		gamma_optimize = gamma_optimize_cvt;
	else
		gamma_optimize = gamma_optimize_odt;

	double gamma = gamma_optimize;
	// energy history
	std::vector<double> energy_history;
	// refinement switch
	bool is_refine_stage = false;
	int cnt_iter = 0;
	int cnt_refine_iter = 0;
	// repeat
	while (cnt_iter < max_iter)
	{
		// BUILD TRIANGULATION BY REMOVING DUPLICATED / CLOSE POINTS AND INSERTING NEW ONES
		rebuild_triangulation(seeds, merging_length);
		
		//		
		if (is_refine_stage && refine_thres > 0) {
			gamma = gamma_refine;
			std::cout << "iter: " << cnt_iter << "(refinement); num of vertices: " << rt_->number_of_vertices() << std::endl;
			//system("pause");
			std::cout << std::endl;
		}
		else {
			gamma = gamma_optimize;
			std::cout << "iter: " << cnt_iter << "; num of vertices: " << rt_->number_of_vertices() << std::endl;
		}
		std::cout << "gamma = " << gamma << std::endl;

		// save mesh
		std::stringstream ss;
		ss << cnt_iter;
		//++cnt_iter;
		std::string fname = "rt_" + ss.str();
		fname += ".mesh";
		write_updated_triangulation(fname);

		// update vertices
		double sum_energy = 0.0;
		double move_sum = 0.0;
		for (auto vit = rt_->finite_vertices_begin(); vit != rt_->finite_vertices_end(); ++vit)
		{
			// skip boundary vertices
			bool is_boundary = false;
			auto vv_cir = rt_->incident_vertices(vit);
			auto vv_end = vv_cir;
			do
			{
				if (rt_->is_infinite(vv_cir))
				{
					is_boundary = true;
					break;
				}
			} while (++vv_cir != vv_end);
			if (is_boundary)
				continue;
			move_sum += update_vertex_i(vit, init_stepsize, alpha, beta, gamma, use_cvt);
			// this one use CVT
			//move_sum += update_vertex_i_alternating(vit, init_stepsize, alpha, beta, gamma, use_cvt);
			Eigen::Vector2d pos(vit->point().point().x(), vit->point().point().y());
			sum_energy += evaluate_F(pos, vit, alpha, beta, gamma, use_cvt);
		}
		sum_energy /= double(rt_->number_of_vertices());
		std::cout << "total movement: " << move_sum << "; sum_of_energy: " << sum_energy << std::endl;
		energy_history.push_back(sum_energy);
		
		// refinement switch
		//if (move_sum < refine_thres)
		//{
		//	is_refine_stage = true;
		//	cnt_refine_iter = 0;
		//}
		if (is_refine_stage)
			cnt_refine_iter++;
		
		cnt_iter++;

		if (cnt_refine_iter > max_refine_iter)
		{
			is_refine_stage = false;
			cnt_refine_iter = 0;
		}
		else if (cnt_iter % 10 == 0) {
			is_refine_stage = true;
		}

		//// measure energy
		//for (auto vit = rt_->finite_vertices_begin(); vit != rt_->finite_vertices_end(); ++vit)
		//{
		//	// skip boundary vertices
		//	bool is_boundary = false;
		//	auto vv_cir = rt_->incident_vertices(vit);
		//	auto vv_end = vv_cir;
		//	do
		//	{
		//		if (rt_->is_infinite(vv_cir))
		//		{
		//			is_boundary = true;
		//			break;
		//		}
		//	} while (++vv_cir != vv_end);
		//	if (is_boundary)
		//		continue;
		//	Eigen::Vector2d pos(vit->point().point().x(), vit->point().point().y());
		//	sum_energy += evaluate_F(pos, vit, alpha, beta, gamma);
		//}
		//sum_energy /= double(rt_->number_of_vertices());
		//std::cout << "total movement: " << move_sum << "; sum_of_energy: " << sum_energy << std::endl;
		//energy_history.push_back(sum_energy);

		// update seeds
		seeds.clear();
		for (auto vit = rt_->finite_vertices_begin(); vit != rt_->finite_vertices_end(); ++vit)
		{
			Point_2 p = vit->point().point();
			double x, y;
			x = p.x(); 
			y = p.y();
			if (p.x() < 0 - offset)
				x = p.x();
			if (p.y() < 0 - offset)
				y = p.y();
			if (p.x() > cols_ + offset)
				x = cols_;
			if (p.y() > rows_ + offset)
				y = rows_;
			seeds.push_back(Point_2(x, y));
		}
	}

	/* OUTPUT */
	write_energy_history("energy_history.list", energy_history);
	write_updated_voronoi("voronoi.mesh");
	write_weighted_vertices("weighted_veritces.array");
}

void CMeshOptimizer::run_optimize_single_pt_as_test()
{
	std::vector<WPoint> seeds;
	for (int ix = 0; ix < x_grid_.size(); ix += 40)
	{
		for (int iy = 0; iy < y_grid_.size(); iy += 40)
		{
			seeds.push_back(WPoint(Point_2(x_grid_[ix], y_grid_[iy]), 0.0));
		}
	}

	//
	rt_->clear();
	rt_->insert(seeds.begin(), seeds.end());

	Point_2 query_pt(80, 110);
	auto vh = rt_->nearest_power_vertex(query_pt);
	Regular_triangulation::Vertex_handle query_vh;
	if ((vh->point().point() - query_pt).squared_length() < 0.0001)
		query_vh = vh;
	else {
		query_vh = rt_->insert(WPoint(query_pt));
	}

	//
	double alpha = 0.001;
	double beta = 0.001;
	double gamma = 1;
	double init_stepsize = 5;

	int cnt_iter = 0;
	int max_iter = 15;
	std::vector<Point_2> trace;
	while (cnt_iter < max_iter)
	{
		// save mesh
		std::stringstream ss;
		ss << cnt_iter;
		++cnt_iter;
		std::string fname = "rt_" + ss.str();
		fname += ".mesh";
		write_updated_triangulation(fname);

		std::cout << "iter: " << cnt_iter << std::endl;
		trace.push_back(query_vh->point().point());
		query_vh->set_new_coord(trace.back());

		if (trace.back().x() < -10 || trace.back().x() > cols_ || trace.back().y() < -10 || trace.back().y() > rows_)
		{
			std::cout << "the pt is moving outside the bounding box..." << std::endl;
			break;
		}
		update_vertex_i(query_vh, init_stepsize, alpha, beta, gamma, false); // do not use cvt for now
	}

	write_single_pt_trace("trace.array", trace);
}

bool CMeshOptimizer::write_input_fields_for_check()
{
	std::cout << "write solid field" << std::endl;
	std::ofstream out_sz_file;
	out_sz_file.open("results\\solid_field.array");
	if (out_sz_file.is_open())
	{
		out_sz_file << z_ << std::endl;
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
		out_sgx_file << zgx_ << std::endl;
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
		out_sgy_file << zgy_ << std::endl;
		out_sgy_file.close();
	}
	else {
		std::cout << "wrong at writing solid grad Y field" << std::endl;
		return false;
	}
	return true;
}

bool CMeshOptimizer::write_updated_triangulation(std::string fname)
{
	std::cout << "write updated triangulation" << std::endl;
	std::ofstream out_mesh_file;
	std::string foldername = "results\\";
	out_mesh_file.open(foldername+fname);
	if (out_mesh_file.is_open())
	{
		//std::map<Regular_triangulation::Vertex_handle, int> vh_id_map;
		//int cnt = 0;
		//for (auto vit = rt_->finite_vertices_begin(); vit != rt_->finite_vertices_end(); ++vit, ++cnt)
		//{
		//	vh_id_map.insert(std::make_pair(vit, cnt));
		//	out_mesh_file << cnt << " " << vit->point().point().x() << " " << vit->point().point().y() << std::endl;
		//}
		//cnt = 0;
		//for (auto fit = rt_->finite_faces_begin(); fit != rt_->finite_faces_end(); ++fit, ++cnt)
		//{
		//	int id0 = vh_id_map[fit->vertex(0)];
		//	int id1 = vh_id_map[fit->vertex(1)];
		//	int id2 = vh_id_map[fit->vertex(2)];
		//	out_mesh_file << cnt << " " << id0 << " " << id1 << " " << id2 << std::endl;
		//}
		//out_mesh_file.close();

		for (auto fit = rt_->finite_faces_begin(); fit != rt_->finite_faces_end(); ++fit)
		{
			out_mesh_file 
				<< fit->vertex(0)->point().point() << " " 
				<< fit->vertex(1)->point().point() << " "
				<< fit->vertex(2)->point().point() << " "
				<< std::endl; 
		}
		out_mesh_file.close();
	}
	else {
		std::cout << "wrong at writing updated triangulation" << std::endl;
		return false;
	}

	return true;
}

bool CMeshOptimizer::write_weighted_vertices(std::string fname)
{
	std::cout << "write weighted vertices" << std::endl;
	std::ofstream out_mesh_file;
	std::string foldername = "results\\";
	out_mesh_file.open(foldername + fname);
	if (out_mesh_file.is_open())
	{
		for (auto vit = rt_->finite_vertices_begin(); vit != rt_->finite_vertices_end(); ++vit)
		{
			//
			bool is_boundary = false;
			auto vf_cir = vit->incident_faces();
			const auto vf_end = vit->incident_faces();
			std::vector<Point_2> omega_polygon;
			do {
				if (rt_->is_infinite(vf_cir))
				{
					is_boundary = true;
					break;
				}

				Point_2 p0 = vf_cir->vertex(0)->point().point();
				Point_2 p1 = vf_cir->vertex(1)->point().point();
				Point_2 p2 = vf_cir->vertex(2)->point().point();
				Point_2 circenter = CGAL::circumcenter(p0, p1, p2);
				omega_polygon.push_back(circenter);
				if (CGAL::area(p0, p1, p2) < 0.000000001)
				{
					std::cerr << "triangle's vertices are almost collinear with tria_area = " << CGAL::area(p0, p1, p2) << std::endl;
					std::cout << p0 << " " << p1 << " " << p2 << std::endl;
					write_updated_triangulation("tangling_mesh.mesh");
				}
			} while (++vf_cir != vf_end);
			if (is_boundary)
				continue;

			omega_polygon.push_back(omega_polygon.front());
			double omega_area = 0.0;
			for (int i = 0; i + 1 < omega_polygon.size(); i++)
			{
				Point_2 p0 = vit->point().point();
				Point_2 p1 = omega_polygon[i];
				Point_2 p2 = omega_polygon[i + 1];
				Point_2 bt = CGAL::barycenter(p0, 1, p1, 1, p2, 1);
				double tria_area = CGAL::area(p0, p1, p2);
				omega_area += tria_area;
			}

			out_mesh_file
				<< vit->point().point() << " "
				<< omega_area
				<< std::endl;
		}
		out_mesh_file.close();
	}
	else {
		std::cout << "wrong at writing updated triangulation" << std::endl;
		return false;
	}

	return true;
}

bool CMeshOptimizer::write_updated_voronoi(std::string fname)
{
	std::cout << "write updated triangulation" << std::endl;
	std::ofstream out_mesh_file;
	std::string foldername = "results\\";
	out_mesh_file.open(foldername + fname);
	if (out_mesh_file.is_open())
	{
		for (auto vit = rt_->finite_vertices_begin(); vit != rt_->finite_vertices_end(); ++vit)
		{
			bool is_boundary = false;

			std::vector<Point_2> cc_list;
			auto vf_cir = vit->incident_faces();
			const auto vf_end = vf_cir;
			do {
				if (rt_->is_infinite(vf_cir))
				{
					is_boundary = true;
					break;
				}
				cc_list.push_back(rt_->weighted_circumcenter(vf_cir));
			} while (++vf_cir != vf_end);
			if (is_boundary)
				continue;

			for (auto cc : cc_list)
			{
				out_mesh_file << cc << " ";
			}
			out_mesh_file << std::endl;
		}
		out_mesh_file.close();
	}
	else {
		std::cout << "wrong at writing updated voronoi" << std::endl;
		return false;
	}

	return true;
}

bool CMeshOptimizer::write_single_pt_trace(std::string fname, std::vector<Point_2> trace)
{
	std::cout << "write single point's trace" << std::endl;
	std::ofstream out_trace_file;
	std::string foldername = "results\\";
	out_trace_file.open(foldername+fname);
	if (out_trace_file.is_open())
	{
		for (auto fit = trace.begin(); fit != trace.end(); ++fit)
		{
			out_trace_file << *fit << std::endl;
		}
		out_trace_file.close();
	}
	else {
		std::cout << "wrong at writing single point's trace" << std::endl;
		return false;
	}

	return true;
}

bool CMeshOptimizer::write_energy_history(std::string fname, std::vector<double> energy_history)
{
	std::cout << "write energy history" << std::endl;
	std::ofstream out_eh_file;
	std::string foldername = "results\\";
	out_eh_file.open(foldername + fname);
	if (out_eh_file.is_open())
	{
		for (auto it = energy_history.begin(); it != energy_history.end(); ++it)
		{
			out_eh_file << *it << std::endl;
		}
		out_eh_file.close();
	}
	else {
		std::cout << "wrong at writing energy history" << std::endl;
		return false;
	}

	return true;
}

void CMeshOptimizer::rebuild_triangulation(std::vector<Point_2> pts, double merging_length)
{
	rt_->clear();
	////// easiest way; but we gonna remove wpts in proximity
	////rt_->insert(wpts.begin(), wpts.end());

	// REMOVE
	std::vector<WPoint> compact_seeds;
	for (auto pt : pts)
	{
		bool has_vertex_in_proximity = false;
		double hpt = dist_transform_field_S(Eigen::Vector2d(pt.x(), pt.y()));
		WPoint wpt = WPoint(pt, 0);
		for (auto seed: compact_seeds)
		{
			double dist = sqrt((wpt.point() - seed.point()).squared_length());
			if (dist < merging_length)
			{
				double hs = dist_transform_field_S(Eigen::Vector2d(seed.point().x(), seed.point().y()));
				if (hs < hpt)
				{
					seed = wpt;
					has_vertex_in_proximity = true;
					break;
				}
			}
		}
		if (!has_vertex_in_proximity) {
			compact_seeds.push_back(wpt);
		}
	}
	rt_->insert(compact_seeds.begin(), compact_seeds.end());

	// INSERT
	insert_seeds_at_edges();
	insert_seeds_at_faces();

}

double CMeshOptimizer::update_vertex_i(Regular_triangulation::Vertex_handle vh, double init_stepsize,
	const double alpha, const double beta, const double gamma, bool use_cvt)
{
	Eigen::Vector2d x_i(vh->point().point().x(), vh->point().point().y());

	/* COMPUTE GRADIENT */
	Eigen::Vector2d grad_F = evaluate_G(x_i, vh, alpha, beta, gamma, use_cvt);

	if (grad_F.norm() < 0.000001)
	{
		std::cerr << "grad is very small" << std::endl;
		return 0.0;
	}
	else 
	{
		/* FEASIBLE REGION */
		double ss = find_feasible_step(x_i, vh, grad_F);

		/* LINE SEARCH */
		double s = init_stepsize < ss ? init_stepsize : ss;
		double s_min = find_steplength_with_wolfe_cond(s, x_i, grad_F, vh, alpha, beta, gamma, use_cvt);

		//double Fcur = evaluate_F(x_i, vh, alpha, beta, gamma);
		//double Fval_min = evaluate_F(x_i + s_min*grad_F, vh, alpha, beta, gamma);
		//std::cout << "Fcur = " << Fcur << ", Fval = " << Fval_min 
		//	<< ", s = " << s_min << ", grad_F = " << grad_F << std::endl;

		//std::cout << "Pos original = " << x_i.transpose() << ", Pos new = " << (x_i + s_min*grad_F).transpose() << std::endl;
		//std::cout << std::endl;

		/* UPDATE X'S POSTION */
		x_i += s_min*grad_F;
		//vh->set_new_coord(Point_2(x_i[0], x_i[1]));
		vh->set_point(WPoint(x_i[0], x_i[1])); //  this will lead to tangling elements
		return s_min*grad_F.norm();
	}
}

double CMeshOptimizer::update_vertex_i_alternating(Regular_triangulation::Vertex_handle vh,
	double init_stepsize, const double alpha, const double beta, const double gamma, bool use_cvt)
{
	Eigen::Vector2d x_i(vh->point().point().x(), vh->point().point().y());

	/* COMPUTE GRADIENT */
	// Do tessellation optimization
	Eigen::Vector2d grad_F_tess = evaluate_G_CVT(x_i, vh);
	//std::cout << x_i.transpose() << ",  " << grad_F_tess.transpose() << std::endl;
	x_i += grad_F_tess;

	//if (grad_F_tess.norm() < 0.000001)
	//{
	//	std::cerr << "grad is very small; we are done at this vertex" << std::endl;
	//	x_i += grad_F_tess;
	//}
	//else
	//{
	//	/* FEASIBLE REGION */
	//	double ss = find_feasible_step(x_i, vh, grad_F_tess);
	//	
	//	/* LINE SEARCH */
	//	double s = init_stepsize < ss ? init_stepsize : ss;
	//	//double s_min = find_steplength_with_wolfe_cond(s, x_i, grad_F_tess, vh, alpha, beta, gamma, use_cvt);
	//	// CONSTANT 
	//	const double c1 = 0.000001;
	//	double Fcur = evaluate_F_constraint(x_i, vh, alpha, beta);
	//	double s_min = 0;
	//	while (s > 0.00001) {
	//		auto new_pos = x_i + s*grad_F_tess;
	//		double Fval = evaluate_F_constraint(new_pos, vh, alpha, beta);
	//		if (Fval <= Fcur + c1*s*grad_F_tess.transpose()*grad_F_tess)
	//		{
	//			s_min = s;
	//			break;
	//		}
	//		s = 0.5*s;
	//	}
	//	/* UPDATE X'S POSTION */
	//	x_i += s_min*grad_F_tess;
	//}

	// Pull back to feasible region
	Eigen::Vector2d grad_F_feasible = evaluate_G_constraint(x_i, vh, alpha, beta);
	if (grad_F_feasible.norm() < 0.000001)
	{
		std::cerr << "grad_F_feasible is very small; we are done at this vertex" << std::endl;
		vh->set_point(WPoint(x_i[0], x_i[1])); 
		return 0.0;
	}
	else
	{
		/* FEASIBLE REGION */
		double ss = find_feasible_step(x_i, vh, grad_F_tess);
		
		/* LINE SEARCH */
		double s = init_stepsize < ss ? init_stepsize : ss;
		
		//double s_min = find_steplength_with_wolfe_cond(s, x_i, grad_F_tess, vh, alpha, beta, gamma, use_cvt);
		// CONSTANT 
		const double c1 = 0.000001;
		double Fcur = evaluate_F_constraint(x_i, vh, alpha, beta);
		double s_min = 0;
		while (s > 0.00001) {
			auto new_pos = x_i + s*grad_F_tess;
			double Fval = evaluate_F_constraint(new_pos, vh, alpha, beta);
			if (Fval <= Fcur + c1*s*grad_F_tess.transpose()*grad_F_tess)
			{
				s_min = s;
				break;
			}
			s = 0.5*s;
		}
		
		/* UPDATE X'S POSTION */
		x_i += s_min*grad_F_tess;
		vh->set_point(WPoint(x_i[0], x_i[1])); 
		return s_min*grad_F_tess.norm();
	}
}

double CMeshOptimizer::evaluate_F(const Eigen::Vector2d & pos, Regular_triangulation::Vertex_handle vh, 
	const double alpha, const double beta, const double gamma, bool use_cvt)
{
	if (use_cvt)
	{
		double f_cvt= evaluate_F_CVT(pos, vh);
		double c = evaluate_F_constraint(pos, vh, alpha, beta);
		return gamma*f_cvt + c;
	}
	else
		return evaluate_F_ODT(pos, vh, alpha, beta, gamma);
}

double CMeshOptimizer::evaluate_F_ODT(const Eigen::Vector2d &pos, Regular_triangulation::Vertex_handle vh,
	const double alpha, const double beta, const double gamma)
{
	Eigen::Matrix2d Rccw90;
	Rccw90 << 0.0, -1.0,
		1.0, 0.0;
	double F = 0.0;
	double E = 0.0; 
	double Ds = 0.0;
	double Dv = 0.0;

	//
	auto vf_cir = vh->incident_faces();
	const auto vf_end = vh->incident_faces();
	double cell_area = 0.0;

	do
	{
		Regular_triangulation::Vertex_handle v0 = vf_cir->vertex(0);
		Regular_triangulation::Vertex_handle v1 = vf_cir->vertex(1);
		Regular_triangulation::Vertex_handle v2 = vf_cir->vertex(2);

		Point_2 pi, pj, pk;
		pi = Point_2(pos[0], pos[1]);
		if (v0 == vh)
		{
			pj = v1->point().point();
			pk = v2->point().point();
		}
		else if (v1 == vh)
		{
			pj = v2->point().point();
			pk = v0->point().point();
		}
		else
		{
			pj = v0->point().point();
			pk = v1->point().point();
		}
		double tria_area = CGAL::area(pi, pj, pk);
		cell_area += tria_area;
		
		// circumcenter and barycenter
		Vector_2 vij_ccw_perp(-(pj - pi).y(), (pj - pi).x());
		Vector_2 vki_ccw_perp(-(pi - pk).y(), (pi - pk).x());
		Point_2 c = pi + (vij_ccw_perp.squared_length()*vki_ccw_perp + vki_ccw_perp.squared_length()*vij_ccw_perp) / (4 * tria_area);
		Point_2 b = Point_2(0,0)+(Vector_2(pi.x(), pi.y()) + Vector_2(pj.x(), pj.y()) + Vector_2(pk.x(), pk.y())) / 3.0;

		Eigen::Vector2d cv(c.x(), c.y());
		Eigen::Vector2d bv(b.x(), b.y());

		// 
		double R2 = (pi - c).squared_length();
		E += tria_area * (R2 + cv.transpose() * (2 * bv - cv)); // E(T_j)
		Dv += tria_area*dist_transform_field_V(Eigen::Vector2d(c.x(), c.y())); // weighing this term by triangle area

	} while (++vf_cir != vf_end);
	
	E /= rows_*cols_; // eliminate the size effect introduced by the coordinate system
	Ds = cell_area*dist_transform_field_S(pos); // weighing this term by cell area
	//Ds = dist_transform_field_S(pos);

	//std::cout << "E = " << E << ", Ds = " << Ds << ", Dv = " << Dv << std::endl;

	return gamma * E+ alpha * Ds+ beta * Dv;
}

double CMeshOptimizer::evaluate_F_CVT(const Eigen::Vector2d & pos, Regular_triangulation::Vertex_handle vh)
{
	// 1 CVT gradient

	// constant
	Eigen::Matrix2d Rccw90;
	Rccw90 << 0.0, -1.0,
		1.0, 0.0;

	Eigen::Vector2d x_i(vh->point().point().x(), vh->point().point().y());

	//
	auto vf_cir = vh->incident_faces();
	const auto vf_end = vh->incident_faces();
	Eigen::Vector2d c_i(0, 0);
	double cell_area = 0.0;
	do {
		Point_2 p0 = vf_cir->vertex(0)->point().point();
		Point_2 p1 = vf_cir->vertex(1)->point().point();
		Point_2 p2 = vf_cir->vertex(2)->point().point();
		double tria_area = CGAL::area(p0, p1, p2);
		Point_2 bt = Point_2(0, 0) + ((p0 - Point_2(0, 0)) + (p1 - Point_2(0, 0)) + (p2 - Point_2(0, 0))) / 3.0;
		cell_area += tria_area;
		c_i += tria_area*Eigen::Vector2d(bt.x(), bt.y());
		if (tria_area < 0.000000001)
		{
			std::cerr << "triangle's vertices are almost collinear with tria_area = " << tria_area << std::endl;
			std::cout << p0 << " " << p1 << " " << p2 << std::endl;
			write_updated_triangulation("tangling_mesh.mesh");
		}
	} while (++vf_cir != vf_end);

	return x_i.transpose()*(cell_area*x_i - c_i);
}

double CMeshOptimizer::evaluate_F_constraint(const Eigen::Vector2d & pos, Regular_triangulation::Vertex_handle vh, const double alpha, const double beta)
{
	// 
	Eigen::Matrix2d Rccw90;
	Rccw90 << 0.0, -1.0,
		1.0, 0.0;
	double Ds = 0.0;
	double Dv = 0.0;

	//
	auto vf_cir = vh->incident_faces();
	const auto vf_end = vh->incident_faces();
	double cell_area = 0.0;
	do
	{
		Regular_triangulation::Vertex_handle v0 = vf_cir->vertex(0);
		Regular_triangulation::Vertex_handle v1 = vf_cir->vertex(1);
		Regular_triangulation::Vertex_handle v2 = vf_cir->vertex(2);

		Point_2 pi, pj, pk;
		pi = Point_2(pos[0], pos[1]);
		if (v0 == vh)
		{
			pj = v1->point().point();
			pk = v2->point().point();
		}
		else if (v1 == vh)
		{
			pj = v2->point().point();
			pk = v0->point().point();
		}
		else
		{
			pj = v0->point().point();
			pk = v1->point().point();
		}
		double tria_area = CGAL::area(pi, pj, pk);
		cell_area += tria_area;

		// circumcenter and barycenter
		Vector_2 vij_ccw_perp(-(pj - pi).y(), (pj - pi).x());
		Vector_2 vki_ccw_perp(-(pi - pk).y(), (pi - pk).x());
		Point_2 c = pi + (vij_ccw_perp.squared_length()*vki_ccw_perp + vki_ccw_perp.squared_length()*vij_ccw_perp) / (4 * tria_area);
		Dv += tria_area*dist_transform_field_V(Eigen::Vector2d(c.x(), c.y())); // weighing this term by triangle area

	} while (++vf_cir != vf_end);

	Ds = cell_area*dist_transform_field_S(pos); // weighing this term by cell area

	return alpha * Ds + beta * Dv;
}

Eigen::Vector2d CMeshOptimizer::evaluate_G(
	const Eigen::Vector2d & pos, Regular_triangulation::Vertex_handle vh, 
	const double alpha, const double beta, const double gamma, bool use_cvt)
{
	if (use_cvt)
	{
		Eigen::Vector2d grad_E = evaluate_G_CVT(pos, vh);
		Eigen::Vector2d grad_C = evaluate_G_constraint(pos, vh, alpha, beta);
		return gamma*grad_E + grad_C;
	}
	else
		return evaluate_G_ODT(pos, vh, alpha, beta, gamma);
}

Eigen::Vector2d CMeshOptimizer::evaluate_G_constraint(const Eigen::Vector2d & pos, Regular_triangulation::Vertex_handle vh, const double alpha, const double beta)
{
	// constant
	Eigen::Matrix2d Rccw90;
	Rccw90 << 0.0, -1.0,
		1.0, 0.0;

	// init gradient
	Eigen::Vector2d grad_Ds;
	grad_Ds.setZero();
	Eigen::Vector2d grad_DvDx;
	grad_DvDx.setZero();

	Eigen::Vector2d x_i(vh->point().point().x(), vh->point().point().y());

	// 2
	grad_Ds = dist_transform_grad_S(x_i);
	assert(grad_Ds.size() == 2);

	//
	auto vf_cir = vh->incident_faces();
	const auto vf_end = vh->incident_faces();
	do {
		Regular_triangulation::Vertex_handle v0 = vf_cir->vertex(0);
		Regular_triangulation::Vertex_handle v1 = vf_cir->vertex(1);
		Regular_triangulation::Vertex_handle v2 = vf_cir->vertex(2);

		Point_2 p0 = v0->point().point();
		Point_2 p1 = v1->point().point();
		Point_2 p2 = v2->point().point();
		double tria_area = CGAL::area(p0, p1, p2);

		if (tria_area < 0.000000001)
		{
			std::cerr << "triangle's vertices are almost collinear with tria_area = " << tria_area << std::endl;
			std::cout << p0 << " " << p1 << " " << p2 << std::endl;
			write_updated_triangulation("tangling_mesh.mesh");
		}

		Eigen::Vector2d x_j, x_k;
		int idx = -1;
		if (v0 == vh)
		{
			idx = 0;
			x_j = Eigen::Vector2d(v1->point().point().x(), v1->point().point().y()) - x_i;
			x_k = Eigen::Vector2d(v2->point().point().x(), v2->point().point().y()) - x_i;
		}
		else if (v1 == vh)
		{
			idx = 1;
			x_j = Eigen::Vector2d(v2->point().point().x(), v2->point().point().y()) - x_i;
			x_k = Eigen::Vector2d(v0->point().point().x(), v0->point().point().y()) - x_i;
		}
		else
		{
			idx = 2;
			x_j = Eigen::Vector2d(v0->point().point().x(), v0->point().point().y()) - x_i;
			x_k = Eigen::Vector2d(v1->point().point().x(), v1->point().point().y()) - x_i;
		}

		// 3
		Eigen::Vector2d grad_DvDx_comp;
		grad_DvDx_comp.setZero();

		// jacobian matrix
		Eigen::Matrix2d Jv, I;
		Jv.setIdentity();
		I.setIdentity();
		Jv += ((x_j - x_i)* (Rccw90*(x_j - x_i)).transpose()
			+ (x_k - x_i)* (Rccw90*(x_k - x_i)).transpose()) / (2 * tria_area);
		Jv += ((x_i - x_j).squaredNorm() * Rccw90*(x_k - x_j) * (Rccw90*(x_j - x_i)).transpose()
			+ (x_i - x_k).squaredNorm() * Rccw90*(x_k - x_j) * (Rccw90*(x_k - x_i)).transpose()) / (8 * tria_area*tria_area);
		Jv += ((x_i - x_j).squaredNorm() + (x_i - x_k).squaredNorm())*I / (4 * tria_area);

		// 
		auto vp_i = rt_->dual(vf_cir);
		Eigen::Vector2d v_i(vp_i.x(), vp_i.y());
		grad_DvDx_comp = Jv*dist_transform_grad_V(v_i);
		grad_DvDx += grad_DvDx_comp;
	} while (++vf_cir != vf_end);

	assert(grad_DvDx.size() == 2);
	assert(grad_Ds.size() == 2);

	grad_DvDx /= double(rows_*cols_);

	Eigen::Vector2d grad_F;
	grad_F = -(alpha*grad_Ds + beta*grad_DvDx);

	//std::cout
	//	<< "grad_E = " << grad_E.transpose() << "\n"
	//	<< "grad_Ds = " << grad_Ds.transpose() << " "
	//	<< " (alpha = " << alpha << ")\n"
	//	<< "grad_DvDx = " << grad_DvDx.transpose() << " "
	//	<< " (beta = " << beta << ")\n";
	//std::cout << "grad_F = " << grad_F.transpose() << " "
	//	<< std::endl;

	return grad_F;
}

Eigen::Vector2d CMeshOptimizer::evaluate_G_ODT(
	const Eigen::Vector2d & pos, 
	Regular_triangulation::Vertex_handle vh, 
	const double alpha, const double beta, const double gamma)
{
	// constant
	Eigen::Matrix2d Rccw90;
	Rccw90 << 0.0, -1.0,
		1.0, 0.0;

	// init gradient
	Eigen::Vector2d grad_Ds;
	grad_Ds.setZero();
	Eigen::Vector2d grad_E;
	grad_E.setZero();
	Eigen::Vector2d grad_DvDx;
	grad_DvDx.setZero();

	Eigen::Vector2d x_i(vh->point().point().x(), vh->point().point().y());

	// 2
	grad_Ds = dist_transform_grad_S(x_i);
	assert(grad_Ds.size() == 2);

	//
	auto vf_cir = vh->incident_faces();
	const auto vf_end = vh->incident_faces();
	double cell_area = 0.0;
	do {
		Regular_triangulation::Vertex_handle v0 = vf_cir->vertex(0);
		Regular_triangulation::Vertex_handle v1 = vf_cir->vertex(1);
		Regular_triangulation::Vertex_handle v2 = vf_cir->vertex(2);

		Point_2 p0 = v0->point().point();
		Point_2 p1 = v1->point().point();
		Point_2 p2 = v2->point().point();
		double tria_area = CGAL::area(p0, p1, p2);

		if (tria_area < 0.000000001)
		{
			std::cerr << "triangle's vertices are almost collinear with tria_area = " << tria_area << std::endl;
			std::cout << p0 << " " << p1 << " " << p2 << std::endl;
			write_updated_triangulation("tangling_mesh.mesh");
		}

		int idx = -1;
		Eigen::Vector2d x_j, x_k;
		if (v0 == vh)
		{
			idx = 0;
			x_j = Eigen::Vector2d(v1->point().point().x(), v1->point().point().y()) - x_i;
			x_k = Eigen::Vector2d(v2->point().point().x(), v2->point().point().y()) - x_i;
		}
		else if (v1 == vh)
		{
			idx = 1;
			x_j = Eigen::Vector2d(v2->point().point().x(), v2->point().point().y()) - x_i;
			x_k = Eigen::Vector2d(v0->point().point().x(), v0->point().point().y()) - x_i;
		}
		else
		{
			idx = 2;
			x_j = Eigen::Vector2d(v0->point().point().x(), v0->point().point().y()) - x_i;
			x_k = Eigen::Vector2d(v1->point().point().x(), v1->point().point().y()) - x_i;
		}

		//// 1
		//grad_E += 1.0 / 3.0 *
		//	((0.5* Rccw90*(x_k - x_j) * (x_k.squaredNorm() + x_j.squaredNorm()))
		//		+
		//		2 * tria_area*x_i);
		grad_E += Rccw90*(x_k - x_j) * (x_k.squaredNorm() + x_j.squaredNorm()) / 6.0;

		//// from Chen ZG's 2014 paper 
		//auto cp_i = rt_->dual(vf_cir);
		//grad_E = -tria_area*Eigen::Vector2d(cp_i.x(), cp_i.y());
		cell_area += tria_area;

		// 3
		Eigen::Vector2d grad_DvDx_comp;
		grad_DvDx_comp.setZero();

		// jacobian matrix
		Eigen::Matrix2d Jv, I;
		Jv.setIdentity();
		I.setIdentity();
		Jv += ((x_j - x_i)* (Rccw90*(x_j - x_i)).transpose()
			+ (x_k - x_i)* (Rccw90*(x_k - x_i)).transpose()) / (2 * tria_area);
		Jv += ((x_i - x_j).squaredNorm() * Rccw90*(x_k - x_j) * (Rccw90*(x_j - x_i)).transpose()
			+ (x_i - x_k).squaredNorm() * Rccw90*(x_k - x_j) * (Rccw90*(x_k - x_i)).transpose()) / (8 * tria_area*tria_area);
		Jv += ((x_i - x_j).squaredNorm() + (x_i - x_k).squaredNorm())*I / (4 * tria_area);

		// 
		auto vp_i = rt_->dual(vf_cir);
		Eigen::Vector2d v_i(vp_i.x(), vp_i.y());
		grad_DvDx_comp = Jv*dist_transform_grad_V(v_i);
		grad_DvDx += grad_DvDx_comp;

		//std::cout << "grad_E: \n " << grad_E << "\ngrad_DvDx: \n " << grad_DvDx << std::endl;
	} while (++vf_cir != vf_end);

	//// from Chen ZG's 2014 paper 
	//grad_E = 2.0*cell_area*x_i/3.0 - 2.0*grad_E / 3.0;

	assert(grad_E.size() == 2);
	assert(grad_DvDx.size() == 2);
	assert(grad_Ds.size() == 2);

	grad_E /= double(rows_*cols_);
	grad_DvDx /= double(rows_*cols_);

	Eigen::Vector2d grad_F;
	grad_F = -(gamma*grad_E + alpha*grad_Ds + beta*grad_DvDx);

	//std::cout
	//	<< "grad_E = " << grad_E.transpose() << "\n"
	//	<< "grad_Ds = " << grad_Ds.transpose() << " "
	//	<< " (alpha = " << alpha << ")\n"
	//	<< "grad_DvDx = " << grad_DvDx.transpose() << " "
	//	<< " (beta = " << beta << ")\n";
	//std::cout << "grad_F = " << grad_F.transpose() << " "
	//	<< std::endl;

	return grad_F;
}

Eigen::Vector2d CMeshOptimizer::evaluate_G_CVT(const Eigen::Vector2d & pos, Regular_triangulation::Vertex_handle vh)
{
	// 1 CVT gradient

	// constant
	Eigen::Matrix2d Rccw90;
	Rccw90 << 0.0, -1.0,
		1.0, 0.0;

	// init gradient
	Eigen::Vector2d grad_E;
	grad_E.setZero();

	Eigen::Vector2d x_i(vh->point().point().x(), vh->point().point().y());

	//
	auto vf_cir = vh->incident_faces();
	const auto vf_end = vh->incident_faces();
	std::vector<Point_2> omega_polygon;
	do {
		Point_2 p0 = vf_cir->vertex(0)->point().point();
		Point_2 p1 = vf_cir->vertex(1)->point().point();
		Point_2 p2 = vf_cir->vertex(2)->point().point();
		Point_2 circenter = CGAL::circumcenter(p0, p1, p2);
		omega_polygon.push_back(circenter);

		if (CGAL::area(p0, p1, p2) < 0.000000001)
		{
			std::cerr << "triangle's vertices are almost collinear with tria_area = " << CGAL::area(p0, p1, p2) << std::endl;
			std::cout << p0 << " " << p1 << " " << p2 << std::endl;
			write_updated_triangulation("tangling_mesh.mesh");
		}
	} while (++vf_cir != vf_end);
	omega_polygon.push_back(omega_polygon.front());

	Eigen::Vector2d c_i(0, 0);
	double omega_area = 0.0;
	for (int i = 0; i+1 < omega_polygon.size(); i++)
	{
		Point_2 p0 = vh->point().point();
		Point_2 p1 = omega_polygon[i];
		Point_2 p2 = omega_polygon[i+1];
		Point_2 bt = CGAL::barycenter(p0, 1, p1, 1, p2, 1);
		double tria_area = CGAL::area(p0, p1, p2);
		omega_area += tria_area;
		c_i += tria_area*Eigen::Vector2d(bt.x(), bt.y());
	}

	//grad_E = 2*cell_area*x_i - 2*c_i;
	grad_E = x_i - c_i/ omega_area; // move half way to c_i

	return -grad_E/2.0;
}

double CMeshOptimizer::find_feasible_step(const Eigen::Vector2d & pos,
	Regular_triangulation::Vertex_handle vh, const Eigen::Vector2d & dir)
{
	Point_2 p(pos[0], pos[1]);
	Vector_2 d(dir[0], dir[1]);
	CGAL::Ray_2<K> ray(p, d);
	double min_step = std::numeric_limits<double>::max();

	/* FIND FEASIBLE REGION */
	auto vv_cir = rt_->incident_vertices(vh);
	const auto vv_end = vv_cir;
	do {
		auto p0 = vv_cir->point().point();
		auto p1 = (++vv_cir)->point().point();
		Line_2 line(p0, p1);
		CGAL::Object obj = CGAL::intersection(ray, line);
		const Point_2* s = CGAL::object_cast<Point_2>(&obj);
		if (s)
		{
			double step = sqrt((*s - p).squared_length()/d.squared_length());
			if (step < min_step)
				min_step = step;
		}
	} while (vv_cir != vv_end);
	return 0.95*min_step;
}

double CMeshOptimizer::find_steplength_with_wolfe_cond(double s, 
	const Eigen::Vector2d &xi, 
	const Eigen::Vector2d &gi,
	const Regular_triangulation::Vertex_handle vh, 
	const double alpha, const double beta, const double gamma, bool use_cvt)
{
	// CONSTANT 
	const double c1 = 0.000001;
	const double c2 = 0.1;
	const bool use_curvature_wolfe_condition = false;

	double Fcur = evaluate_F(xi, vh, alpha, beta, gamma, use_cvt);
	double s_min = 0;
	//
	while(s > 0.00001) {
		auto new_pos = xi + s*gi;
		double Fval = evaluate_F(new_pos, vh, alpha, beta, gamma, use_cvt);

		//std::cout << "s = " << s << ", pos = " << new_pos.transpose() << std::endl;
		//std::cout << "Fcur = " << Fcur << ", Fval = " << Fval << std::endl;

		if (Fval <= Fcur + c1*s*gi.transpose()*gi)
		{
			if (!use_curvature_wolfe_condition)
			{
				s_min = s;
				break;
			}
			// Curvature Wolfe condition
			else {
				auto gii = evaluate_G(xi + s_min*gi, vh, alpha, beta, gamma, use_cvt);
				if (-gi.transpose()*gii <= -c2*gi.transpose()*gi)
				{
					s_min = s;
					break;
				}
			}
		}
		s = 0.5*s;
	}

	// check
	CGAL::Polygon_2<K> polygon;
	auto vv_cir = rt_->incident_vertices(vh);
	const auto vv_end = vv_cir;
	do {
		polygon.push_back(vv_cir->point().point());
	} while (++vv_cir != vv_end);

	auto xnew = xi + s_min*gi;
	Point_2 pi(xnew[0], xnew[1]);
	assert(polygon.bounded_side(pi) == CGAL::ON_BOUNDED_SIDE);
	return s_min;
}

void CMeshOptimizer::insert_seeds_at_edges()
{
	std::vector<WPoint> seeds_to_insert;

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

		//auto psrc = Point_2(2.3, 147.1);
		//auto ptgt = Point_2(45.2, 148.0);
		//
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
			seeds_to_insert.push_back(WPoint(p, 0));
		}
	}

	std::cout << "# new inserted seeds: " << seeds_to_insert.size() << std::endl;

	// INSERT
	rt_->insert(seeds_to_insert.begin(), seeds_to_insert.end());

}

void CMeshOptimizer::insert_seeds_at_faces()
{
	//std::vector<WPoint> seeds_to_insert;

	//// for each face in the new built triangulation, 
	//// if the face contains a local minimum, then we shall add a point there to describe it

	//for (auto fit = rt_->finite_faces_begin(); fit != rt_->finite_faces_end(); ++fit)
	//{
	//	auto fh = eit->first;
	//	auto vh = fh->vertex(eit->second);
	//	auto psrc = fh->vertex(fh->ccw(eit->second))->point().point();
	//	auto ptgt = fh->vertex(fh->cw(eit->second))->point().point();

	//	//auto psrc = Point_2(2.3, 147.1);
	//	//auto ptgt = Point_2(45.2, 148.0);
	//	//
	//	// sampling on edge(psrc, ptgt)
	//	Vector_2 vdir = ptgt - psrc;
	//	if (vdir.squared_length() < 0.00001)
	//	{
	//		// should not have these edges because of merging is done ahead of insersion
	//		// skip very short edges
	//		std::cerr << "very short edge" << std::endl;
	//		//continue;
	//	}
	//	// we only sample those in the middle, because the two end should have grad = 0
	//	double step = 0.1;
	//	double min_grad_q_norm = 4;
	//	double min_t = -1;
	//	int cnt_turning_sign = 0;
	//	double sign = -1;
	//	double h0 = dist_transform_field_S(Eigen::Vector2d(psrc.x(), psrc.y()));
	//	for (double t = 0.05; t < 0.95; t += step)
	//	{
	//		Point_2 q = psrc + t*vdir;
	//		double hq = dist_transform_field_S(Eigen::Vector2d(q.x(), q.y()));
	//		if (sign * hq < 0)
	//		{
	//			sign = -sign;
	//			cnt_turning_sign++;
	//		}
	//		Eigen::Vector2d grad_q = dist_transform_grad_S(Eigen::Vector2d(q.x(), q.y()));
	//		if (hq < 0 && grad_q.norm() < min_grad_q_norm && cnt_turning_sign == 2)
	//		{
	//			min_t = t;
	//			min_grad_q_norm = grad_q.norm();
	//		}
	//	}
	//	if (min_t > 0 && cnt_turning_sign > 3)
	//	{
	//		auto p = psrc + min_t*vdir;
	//		double hinsert = dist_transform_field_S(Eigen::Vector2d(p.x(), p.y()));
	//		seeds_to_insert.push_back(WPoint(p, 0));
	//	}
	//}

	//std::cout << "# new inserted seeds: " << seeds_to_insert.size() << std::endl;

	//// INSERT
	//rt_->insert(seeds_to_insert.begin(), seeds_to_insert.end());
}

void CMeshOptimizer::boundary_sampling(Point_2 psrc, Point_2 ptgt, std::vector<Point_2> &out_samples)
{
	// find the negative interval 
	std::vector<std::pair<Point_2, Point_2>> intervals;
	// north boundary for rectangle
	Vector_2 vdir = ptgt - psrc;
	// we only sample those in the middle, because the two end should have grad = 0
	double step = 0.01;
	int cnt_turning_sign = 0;
	double h0 = dist_transform_field_S(Eigen::Vector2d(psrc.x(), psrc.y()));
	double sign = h0/fabs(h0);
	Point_2 start(-1, -1);
	bool positive_start = true;
	if (sign < 0)
	{
		start = psrc;
		positive_start = false;
	}
	for (double t = 0.0; t < 1.0; t += step)
	{
		Point_2 q = psrc + t*vdir;
		double hq = dist_transform_field_S(Eigen::Vector2d(q.x(), q.y()));
		if (sign * hq < 0)
		{
			sign = -sign;
			cnt_turning_sign++;
			// from positive to negative
			if (positive_start) {
				if (cnt_turning_sign % 2 == 1)
				{
					start = q;
				}
				else {
					//switch_off = q;
					//std::cout << "[ " << start << " -> " << q << "], " << std::endl;
					intervals.push_back(std::make_pair(start, q));
					start = Point_2(-1, -1);
				}
			}
			// from negative start to negative
			else {
				if (cnt_turning_sign % 2 == 0)
				{
					start = q;
				}
				else {
					//std::cout << "[ " << start << " -> " << q << "], " << std::endl;
					intervals.push_back(std::make_pair(start, q));
					start = Point_2(-1, -1);
				}
			}
		}
	}

	//std::cout << std::endl;

	// find lowest point in each interval
	std::vector<Point_2> samples;
	for (auto an_interval : intervals)
	{
		Vector_2 vdir = an_interval.second - an_interval.first;
		double step = 0.1;
		double min_t = -1;
		double min_hq = std::numeric_limits<double>::max();
		for (double t = 0.0; t < 1.0; t += step)
		{
			Point_2 q = an_interval.first + t*vdir;
			double hq = dist_transform_field_S(Eigen::Vector2d(q.x(), q.y()));
			if (hq < min_hq)
			{
				min_t = t;
				min_hq = hq;
			}
		}
		samples.push_back(an_interval.first + min_t*vdir);
	}
	out_samples = samples;
}

void CMeshOptimizer::compute_field_gradients(
	const Eigen::MatrixXd &zmap,
	Eigen::MatrixXd &zgradx, 
	Eigen::MatrixXd &zgrady)
{
	zgradx.resize(zmap.rows(), zmap.cols());
	zgradx.setZero();
	
	zgrady.resize(zmap.rows(), zmap.cols());
	zgrady.setZero();

	double gx, gy;
	for (int ix = 1; ix+1 < x_grid_.size(); ix++)
	{
		for (int iy = 1; iy+1 < y_grid_.size(); iy++)
		{
			// central difference
			gy = (zmap(iy + 1, ix - 1) - zmap(iy - 1, ix - 1)
				+ zmap(iy + 1, ix    ) - zmap(iy - 1, ix    )
				+ zmap(iy + 1, ix + 1) - zmap(iy - 1, ix + 1)) / 3.0;

			gx = (zmap(iy - 1, ix + 1) - zmap(iy - 1, ix - 1)
				+ zmap(iy    , ix + 1) - zmap(iy    , ix - 1)
				+ zmap(iy + 1, ix + 1) - zmap(iy + 1, ix - 1)) / 3.0;

			zgradx(iy, ix) = gx;
			zgrady(iy, ix) = gy;
		}
	}
}

double CMeshOptimizer::get_interp_val_from_2dgrid(const Eigen::Vector2d & pos, const Eigen::MatrixXd & map)
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

