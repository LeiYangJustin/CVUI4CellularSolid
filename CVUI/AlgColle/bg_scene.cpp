#include "bg_scene.h"
#include <ctime>
#include <ratio>
#include <chrono>

void BgScene::construct_bg_info(double cols, double rows,
	Eigen::MatrixXd & z, Eigen::MatrixXd & zgx, Eigen::MatrixXd & zgy,
	std::vector<Point_2> in_edgepts)
{
	cols_ = (cols - 1) / 2.0;
	rows_ = (rows - 1) / 2.0;
	
	Vector_2 shift_vec(-(cols - 1) / 2.0, -(rows - 1) / 2.0);
	for (int i = 0; i < cols - 1; i++)
	{
		for (int j = 0; j < rows - 1; j++)
		{
			Point_2 pxy(i, j);
			z_map_.insert(std::make_pair(pxy + shift_vec, z(i, j)));
			zgx_map_.insert(std::make_pair(pxy + shift_vec, zgx(i, j)));
			zgy_map_.insert(std::make_pair(pxy + shift_vec, zgy(i, j)));
			//z_map_.insert(std::make_pair(pxy + shift_vec, z(j, i)));
			//zgx_map_.insert(std::make_pair(pxy + shift_vec, zgx(j, i)));
			//zgy_map_.insert(std::make_pair(pxy + shift_vec, zgy(j, i)));
		}
	}
	for (int i = 0; i <in_edgepts.size(); ++i)
	{
		in_edgepts[i] += shift_vec;
	}

	x_interval_.clear();
	y_interval_.clear();
	for (int i = 0; i < cols - 1; i++)
	{
		x_interval_.push_back(i + shift_vec.x());
	}
	for (int j = 0; j < rows - 1; j++)
	{
		y_interval_.push_back(j + shift_vec.y());
	}

	// build kd-tree for quick query of edge pts
	tree_search_edgepts_tool_.build_tree(in_edgepts);
	has_edgepts_kdtree_ = true;

	//std::vector<Point_2> NNpts;
	//tree_search_edgepts_tool_.get_fr_nn_from_tree(Point_2(0, 0), NNpts, 100);
}

unsigned BgScene::bg_optimize_all_for_extraction(FT wstep, FT xstep,
	unsigned max_newton_iters,
	FT epsilon, unsigned max_iters,
	std::ostream& out)
{
	out << "#vertices: " << m_rt.number_of_vertices() << std::endl;

	// initialize
	FT norm = optimize_positions_via_lloyd(true);
	//FT norm = 1000;
	FT xthreshold = compute_position_threshold(epsilon);
	xthreshold = 1.0;
	out << "Threshold: " << xthreshold << std::endl;

	std::string fname = "synthss_";
	std::string fname_increment = fname + "0";
	std::string affix = ".mesh";
	write_updated_triangulation(fname_increment + affix, m_rt);

	// loop
	unsigned iters = 0;
	while (iters < max_iters)
	{
		std::cout << "bg_optimize_all_for_extraction" << std::endl;
		iters++;
		//norm = bg_optimize_positions_via_combined_gradient_ascent(xstep, true);
		norm = bg_optimize_positions_via_lagrangian_multiplier();

		std::cout << "Iter: " << iters << "; Norm: " << norm << std::endl;
		if (norm <= xthreshold) break;

		// INSERT PTS
		std::vector<Point_2> new_sites;
		std::vector<double> new_weights;
		//insert_seeds_at_faces();
		double merging_length = 6;
		insert_seeds_at_edges(new_sites, merging_length);
		remove_seeds_in_proximity(new_sites, 0.5*merging_length);

		std::cout << "old sites: " << m_rt.number_of_vertices() 
			<< "new sites: " << new_sites.size() 
			<< "\n" << std::endl;

		for (int i = 0; i < new_sites.size(); i++)
			new_weights.push_back(0.0);
		
		// reset
		construct_triangulation(new_sites, new_weights);
		// write new triangulation
		fname_increment = fname + std::to_string(iters);
		write_updated_triangulation(fname_increment + affix, m_rt);
	}

	return -1;
}

FT BgScene::bg_optimize_positions_via_combined_gradient_ascent(FT timestep,
	bool update)
{
	alpha_ = 1;
	solver_type_ = FIELD_SOLVER;

	std::vector<Point_2> points;
	collect_visible_points(points);
	// gradient
	std::vector<Vector_2> bg_gradient;
	bg_compute_position_gradient(bg_gradient, 1.0, get_solvertype());

	if (timestep <= 0.0)
	{
		//double mean_capacity = compute_mean(m_capacities);
		//double max_alpha = 1.0/ mean_capacity;
		double max_alpha = 1.0;
		int max_iter_linesearch = 20;
		BG_LSPositions bg_line_search(this, max_iter_linesearch, max_alpha);

		FT step = bg_line_search.run_bt(points, bg_gradient);
		std::cout << "linesearch: step " << step << std::endl;
		//
		bg_gradient.clear();
		bg_compute_position_gradient(bg_gradient, 1.0, get_solvertype());

		return compute_norm(bg_gradient);
	}

	return compute_norm(bg_gradient);
}

FT BgScene::bg_optimize_positions_via_lagrangian_multiplier()
{
	std::cout << "bg_optimize_positions_via_lagrangian_multiplier" << std::endl;

	double norm = 0.0;
	alpha_ = 0.5;

	// init miu and lambda
	solver_type_ = ALM_SOLVER;
	alm_lambda_.clear();
	for (int i = 0; i < m_vertices.size(); ++i)
	{
		alm_lambda_.push_back(0.5);
	}

	//this->toggle_connectivity();
	std::vector<Vector_2> bg_gradient;
	std::vector<Point_2> points;
	for (int k = 0; k < 5; k++)
	{
		assert(alm_lambda_.size() > 0);
		// increase lambda to enforce the constraint
		for (int i = 0; i < alm_lambda_.size(); ++i)
		{
			alm_lambda_[i] = pow(3, k); // 0.5*3^4 = 28.8325195313
		}

		// gradient descent to search for update
		points.clear();
		collect_visible_points(points);

		bg_gradient.clear();
		bg_compute_position_gradient_alm_solver(bg_gradient, 1.0, alm_lambda_);

		// normalize bg_gradient to make them local
		for (auto & g : bg_gradient)
		{
			g /= sqrt(g.squared_length());
		}



		double max_alpha = 5.0;
		int max_iter_linesearch = 20;
		BG_LSPositions bg_line_search(this, max_iter_linesearch, max_alpha);
		FT step = bg_line_search.run_bt(points, bg_gradient);
		std::cout << "linesearch: step " << step << std::endl;

		// grad at new pos
		bg_gradient.clear();
		bg_compute_position_gradient(bg_gradient, 1.0, get_solvertype());
		std::cout << "grad_norm: " << compute_norm(bg_gradient) << std::endl;
		if (compute_norm(bg_gradient) < 0.000001)
			break;
	}
	
	//bg_gradient.clear();
	//bg_compute_position_gradient(bg_gradient, 1.0, get_solvertype());
	//this->toggle_connectivity();
	return compute_norm(bg_gradient);
}

void BgScene::init_alm_parameters()
{
	// init miu and lambda
	solver_type_ = ALM_SOLVER;
	alm_miu_ = 1.0;
	alm_lambda_;
	for (int i = 0; i < m_vertices.size(); ++i)
	{
		alm_lambda_.push_back(1.0);
	}
}

void BgScene::bg_compute_position_gradient(std::vector<Vector_2>& gradient, const FT coef, SolverType stype)
{
	switch (stype)
	{
	case BgScene::CVT_SOLVER:
		bg_compute_position_gradient_cvt_solver(gradient, coef);
		break;
	case BgScene::FIELD_SOLVER:
		bg_compute_position_gradient_field_solver(gradient, coef);
		break;
	case BgScene::COMBINED_CVT_FIELD_SOLVER:
		bg_compute_position_gradient_combined_cvt_field_solver(gradient, coef);
		break;
	case BgScene::ALM_SOLVER:
		bg_compute_position_gradient_alm_solver(gradient, coef, alm_lambda_);
		break;
	default:
		break;
	}
}

void BgScene::bg_compute_position_gradient_cvt_solver(std::vector<Vector_2>& gradient, const FT coef)
{
	compute_position_gradient(gradient, coef);
	double mean_capacity = compute_mean(m_capacities);
	for (unsigned i = 0; i < gradient.size(); ++i)
	{
		gradient[i] /= mean_capacity;
	}
}

void BgScene::bg_compute_position_gradient_field_solver(std::vector<Vector_2>& gradient, const FT coef)
{
	gradient.clear();

	// bg constraint
	for (unsigned i = 0; i < m_vertices.size(); ++i)
	{
		Vertex_handle vh = m_vertices[i];
		// Solid phase
		Eigen::Vector2d grad_Ds;
		grad_Ds.setZero();
		// -1.0 indicates negative gradient; 
		grad_Ds = -1.0*dist_transform_grad_S(vh->point().point());
		// coef is used to scale the gradient or to get positive gradient when set to -1.0
		gradient.push_back(Vector_2(coef*grad_Ds[0], coef*grad_Ds[1]));
	}
}

void BgScene::bg_compute_position_gradient_combined_cvt_field_solver(std::vector<Vector_2>& gradient, const FT coef)
{
	gradient.clear();

	std::vector<Vector_2> cvt_gradient, field_gradient;
	bg_compute_position_gradient_cvt_solver(cvt_gradient, coef);
	bg_compute_position_gradient_field_solver(field_gradient, coef);

	assert(alpha_ >= 0);
	assert(m_vertices.size() == cvt_gradient.size());
	assert(m_vertices.size() == field_gradient.size());
	for (int i = 0; i < m_vertices.size(); ++i)
	{
		// alpha_ is a scaling factor; must be positive
		Vector_2 gi = alpha_* cvt_gradient[i] + field_gradient[i];

		//std::cout << "G_cvt = " << alpha_* cvt_gradient[i]
		//	<< "G_f = " << field_gradient[i]
		//	<< std::endl;

		gradient.push_back(gi);
	}
}

void BgScene::bg_compute_position_gradient_alm_solver(
	std::vector<Vector_2>& gradient, const FT coef, const std::vector<double> lambda)
{
	std::vector<Vector_2> cvt_gradient, field_gradient;
	bg_compute_position_gradient_cvt_solver(cvt_gradient, coef);
	bg_compute_position_gradient_field_solver(field_gradient, coef);

	assert(lambda.size() == m_vertices.size());
	assert(alpha_ >= 0);
	for (int i = 0; i < m_vertices.size(); ++i)
	{
		// pgi = -g_cvt - lambda[i]*\nabla(c(x)[i])
		// c(x)[i] = field_gradient
		
		//// feasible moving direction
		auto g = cvt_gradient[i];
		project_gradient_to_feasible_direction(g, m_vertices[i]->get_position());
		Vector_2 pgi = alpha_ * g + 2 * lambda[i] * field_gradient[i];

		//Vector_2 pgi = alpha_ * cvt_gradient[i] + 2*lambda[i]*field_gradient[i];

		//std::cout << "Pt: " << m_vertices[i]->point().point()
		//	<< ", g = " << g 
		//	<< ", G_cvt = " << cvt_gradient[i]
		//	<< ", G_f = " << field_gradient[i]
		//	<< ", pgi = " << pgi 
		//	<< std::endl;

		// this is negative gradient (when coef > 0)
		gradient.push_back(pgi);

		//std::cout << pgi << std::endl;
	}
}

FT BgScene::bg_compute_energy(SolverType stype)
{
	switch (stype)
	{
	case BgScene::CVT_SOLVER:
		return bg_compute_energy_wcvt_solver();
		break;
	case BgScene::FIELD_SOLVER:
		return bg_compute_energy_field_solver();
		break;
	case BgScene::COMBINED_CVT_FIELD_SOLVER:
		return bg_compute_energy_combined_cvt_field_solver();
		break;
	case BgScene::ALM_SOLVER:
		return bg_compute_energy_alm_solver(alm_lambda_);
		break;
	default:
		return bg_compute_energy_wcvt_solver();
		break;
	}
}

FT BgScene::bg_compute_energy_wcvt_solver()
{
	return compute_wcvt_energy();
}

FT BgScene::bg_compute_energy_field_solver()
{
	FT fVal = 0;

	std::vector<Vector_2> field_gradient;
	bg_compute_position_gradient_field_solver(field_gradient, 1.0);

	// bg constraint energy
	for (unsigned i = 0; i < m_vertices.size(); ++i)
	{
		Vertex_handle vh = m_vertices[i];
		Point_2	pi = vh->point().point();

		auto vf_cir = m_rt.incident_faces(vh);
		const auto vf_end = vf_cir;
		double cell_area = vh->compute_area();
		double Ds = 0.0;

		// weighing this term by cell area
		
		Ds = dist_transform_field_S(pi); 
		fVal += Ds;
	}

	return -fVal; // reverse the sign to be consistent with de Goes's code
}

FT BgScene::bg_compute_energy_combined_cvt_field_solver()
{
	FT fVal = 0.0;
	assert(alpha_ >= 0);
	fVal = alpha_ * bg_compute_energy_wcvt_solver() + bg_compute_energy_field_solver();
	return fVal;
}

FT BgScene::bg_compute_energy_alm_solver(const std::vector<double> lambda)
{
	assert(alpha_ >= 0);
	
	FT fVal = 0.0;

	std::vector<Vector_2> field_gradient;
	bg_compute_position_gradient_field_solver(field_gradient, 1.0);

	double cvtE = compute_wcvt_energy();
	fVal = alpha_ * cvtE;
	
	for (int i = 0; i < m_vertices.size(); ++i) {
		fVal -= lambda[i] * field_gradient[i].squared_length();
	}

	// 
	fVal += lambda[0] * bg_compute_energy_field_solver();

	return fVal;
}

//double BgScene::get_interp_val_from_kdtree(const Point_2 &pos, std::map<Point_2, FT> &map, const double coef)
//{
//	if (!has_kdtree)
//		return -1;
//
//	return -1;
//}

void BgScene::transform_by_scaled_eigen_vectors(
	Vector_2 &gradient,	const std::vector<Point_2> datapts, const Point_2 pt_cntr)
{
	//std::cout << "scaled_eigen_vectors_of_pts" << std::endl;
	int nData = datapts.size();

	// center/mean
	Eigen::Vector2d cntr(0, 0);
	if (pt_cntr == Point_2())
	{
		// compute mean
		for (unsigned i = 0; i < nData; ++i)
			cntr += Eigen::Vector2d(datapts[i].x(), datapts[i].y());
		cntr /= double(nData);
	}
	else {
		cntr = Eigen::Vector2d(pt_cntr.x(), pt_cntr.y());
	}

	// weighted covariance
	Eigen::Matrix2d Cmat(2, 2);
	Cmat.setZero();
	for (unsigned i = 0; i < nData; ++i)
	{
		Eigen::Vector2d v = Eigen::Vector2d(datapts[i].x(), datapts[i].y()) - cntr;
		//double weight = exp(4.0 * v.squaredNorm() / 9.0 + 0.000001);
		double weight = 1.0;
		Cmat += weight*(v*v.transpose());
	}

	// eigen_decompose covariance matrix
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(Cmat);
	if (eigensolver.info() != Eigen::Success) {
		std::cerr << "eigen decomposition failed" << std::endl;
	}
	Eigen::Matrix2d scaleMat;
	scaleMat.setZero();
	scaleMat(0, 0) = eigensolver.eigenvalues()[0]/ eigensolver.eigenvalues()[1];
	scaleMat(1, 1) = 1.0;
	Eigen::Matrix2d rotMat;
	rotMat = eigensolver.eigenvectors();
	// transformation
	Eigen::Vector2d g(gradient.x(), gradient.y());
	Eigen::Vector2d gt = rotMat.transpose()*scaleMat*rotMat*g;
	gradient = Vector_2(gt[0], gt[1]);

	if (isnan(gt[0]) || isnan(gt[1]))
	{
		//std::cerr << gradient << std::endl;
		std::cout << "num of NNpts: " << nData << std::endl;
		std::cout << "eigenvectors: \n" << rotMat << std::endl;
		std::cout << "eigenvalues: \n" << scaleMat << std::endl;
		std::cout << "before transform: \n" << g.transpose() << std::endl;
		std::cout << "after transform: \n" << gt.transpose() << std::endl;
	}
}

void BgScene::project_gradient_to_feasible_direction(Vector_2 & gradient, const Point_2 pt)
{
	//std::cout << "project_gradient_to_feasible_direction" << std::endl;

	// find NN pts
	std::vector<Point_2> nnpts;
	double sq_radius = 400.0;
	tree_search_edgepts_tool_.get_fr_nn_from_tree(pt, nnpts, sq_radius);
	if (nnpts.size() > 3) {
		transform_by_scaled_eigen_vectors(gradient, nnpts, pt);
	}
}

double BgScene::get_interp_val_from_2dgrid(const Point_2 &pos, std::map<Point_2, FT> &map, const double coef)
{
	double xl, xu, yl, yu;
	get_interval(pos.x(), xl, xu, x_interval_);
	get_interval(pos.y(), yl, yu, y_interval_);

	Point_2 p00, p01, p11, p10;
	p00 = Point_2(xl, yl);
	p01 = Point_2(xu, yl);
	p10 = Point_2(xu, yu);
	p11 = Point_2(xl, yu);

	//std::cout << pos << ": "<< xl << ", " << xu << "; " << yl << ", " << yu << std::endl;

	// bilinear interpolation
	double s = pos.x() - xl;
	double t = pos.y() - yl;

	// check this: x -> cols, 
	FT v00, v01, v11, v10;
	if (map.find(p00) != map.end())
	{
		v00 = map[p00];
	}
	else {
		std::cout << p00 << std::endl;
	}

	if (map.find(p01) != map.end())
	{
		v01 = map[p01];
	}
	else {
		std::cout << p01 << std::endl;
	}

	if (map.find(p11) != map.end())
	{
		v11 = map[p11];
	}
	else {
		std::cout << p11 << std::endl;
	}

	if (map.find(p10) != map.end())
	{
		v10 = map[p10];
	}
	else {
		std::cout << p10 << std::endl;
	}

	//std::cout << "p00 " << p00 << ": v00 " << v00 << "\n"
	//	<< "p01 " << p01 << ": v01 " << v01 << "\n"
	//	<< "p11 " << p11 << ": v11 " << v11 << "\n"
	//	<< "p10 " << p10 << ": v10 " << v10 << "\n"
	//	<< std::endl;

	double v_interp = (1 - s)*(1 - t)*v00 + (1 - s)*t*v01 + s*t*v11 + s*(1 - t)*v10;

	return coef * v_interp;
}

void BgScene::get_interval(double x, double & xl, double & xu, const std::vector<double> &interval)
{
	int idl, idu;
	idl = -1;
	idu = -1;
	double dl, du;
	dl = std::numeric_limits<double>::max();
	du = std::numeric_limits<double>::max();
	for (int i = 0; i < interval.size(); ++i)
	{
		// find upper
		double d = interval[i] - x;
		if (d < 0 && fabs(d) < dl)
		{
			dl = fabs(d);
			idl = i;
		}
		if (d > 0 && fabs(d) < du)
		{
			du = fabs(d);
			idu = i;
		}
	}

	// in between
	if (idl != -1 && idu != -1)
	{
		xl = interval[idl];
		xu = interval[idu];
		return;
	}
	// no lower
	else if (idl == -1)
	{
		xl = interval[idu];
		xu = interval[idu+1];
	}
	// no upper
	else if (idu == -1)
	{
		xl = interval[idl-1];
		xu = interval[idl];
	}
}


void BgScene::remove_seeds_in_proximity(std::vector<Point_2> &new_sites, double merging_length)
{
	// REMOVE
	std::vector<RT::Point_2> compact_seeds;
	for (auto vit = m_rt.finite_vertices_begin(); vit != m_rt.finite_vertices_end(); ++vit)
	{
		bool has_vertex_in_proximity = false;
		double hpt = dist_transform_field_S(vit->point().point());
		//WPoint wpt = WPoint(pt, 0);
		for (auto seed : compact_seeds)
		{
			double dist = sqrt((vit->point().point() - seed).squared_length());
			if (dist < merging_length)
			{
				double hs = dist_transform_field_S(seed);
				if (hs < hpt)
				{
					seed = vit->point().point();
					has_vertex_in_proximity = true;
					break;
				}
			}
		}
		if (!has_vertex_in_proximity) {
			compact_seeds.push_back(vit->point().point());
		}
	}

	//new_sites.clear();
	//new_sites = compact_seeds;
	new_sites.insert(new_sites.end(), compact_seeds.begin(), compact_seeds.end());
}

void BgScene::insert_seeds_at_edges(std::vector<Point_2> &new_sites, double merging_length)
{
	std::vector<Point_2> seeds_to_insert;
	std::vector<std::pair<Point_2, std::pair<RT::Face_handle, int>>> seed_edge_list_to_insert;
	// for each edge in the new built triangulation, 
	// if the edge passes a local minimum, then we shall add a point there to describe it

	for (auto eit = m_rt.all_edges_begin(); eit != m_rt.all_edges_end(); ++eit)
	{
		// skip infinite edges
		if (m_rt.is_infinite(eit))
			continue;

		auto fh = eit->first;
		auto vh = fh->vertex(eit->second);
		auto psrc = fh->vertex(fh->ccw(eit->second))->point().point();
		auto ptgt = fh->vertex(fh->cw(eit->second))->point().point();

		//psrc = RT::Point(13.3, 32.1);
		//ptgt = RT::Point(0.24, 28.2);

		// sampling on edge(psrc, ptgt)
		Vector_2 vdir = ptgt - psrc;
		if (vdir.squared_length() < 0.00001)
		{
			// should not have these edges because of merging is done ahead of insersion
			// skip very short edges
			std::cerr << "very short edge" << std::endl;
			continue;
		}

		// we only sample those in the middle, because the two end should have grad = 0
		double step = 0.2;
		double min_grad_q_norm = 4;
		double min_t = -1;
		int cnt_turning_sign = 0;
		double sign = -1;
		double h0 = dist_transform_field_S(psrc);
		for (double t = 0.05; t < 0.95; t += step)
		{
			Point_2 q = psrc + t*vdir;
			double hq = dist_transform_field_S(q);
			if (sign * hq < 0)
			{
				sign = -sign;
				cnt_turning_sign++;
			}
			Eigen::Vector2d grad_q = dist_transform_grad_S(q);
			if (hq < 0 && grad_q.norm() < min_grad_q_norm && cnt_turning_sign == 2)
			{
				min_t = t;
				min_grad_q_norm = grad_q.norm();
			}
		}
		double disp = min_t*sqrt(vdir.squared_length());
		if (min_t > 0 && cnt_turning_sign > 3 && disp > merging_length)
		{
			auto p = psrc + min_t*vdir;
			double hinsert = dist_transform_field_S(p);
			std::cout << "inserted: " << p << ", t = " << min_t << ", disp = " << disp << std::endl;
			seeds_to_insert.push_back(p);
			seed_edge_list_to_insert.push_back(std::make_pair(p, std::make_pair(eit->first, eit->second)));
		}
	}
	//std::cout << "# new inserted sites: " << seeds_to_insert.size() << std::endl;
	//new_sites = seeds_to_insert;

	for (auto seed_edge : seed_edge_list_to_insert) {
		auto p = seed_edge.first;
		auto fh = seed_edge.second.first;
		auto eid = seed_edge.second.second;
		m_rt.insert_in_edge(RT::Weighted_point(p, 0), fh, eid);
	}
}

void BgScene::insert_seeds_at_faces()
{
	std::cout << "insert_seeds_at_faces under construction" << std::endl;
}