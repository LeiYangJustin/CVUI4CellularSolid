#include "synthesis.h"

#include<Eigen/SparseQR>

#include "geo_calculator.h"

bool CMeshSynthesizer::synthesize()
{
	// MAX ITER FOR THE GLOBAL LOOP
	const int MAX_ITER = 50;

	// PARAMETERS FOR SYNTHESIS
	int verbose = 0;
	int stepX = 0.0;
	int stepW = 0.0;
	int epsilon = 1.0;
	int frequency = 0;
	int max_newton_iters = 500;
	int max_opt_iters = 500;


	// INIT THE TARGET DOMAIN WITH PTS
	if (!has_example_)
	{
		std::cout << "pls set the exemplar" << std::endl;
		return false;
	}
	if (!has_domain_)
	{
		std::cout << "pls set the target domain to be filled" << std::endl;
		return false;
	}

	// INIT THE SCENE OBJECT
	int nb = 100;
	scene_.generate_random_sites(nb);

	/* need some codes to assign capacity constraints from the example */

	// REPEAT
	int syn_iter = 0;
	while (syn_iter++ < MAX_ITER)
	{
		// SCENE_OPTIMIZE
		scene_.optimize_all(stepW, stepX, max_newton_iters, epsilon, max_opt_iters, std::cout);
		//scene_.get_rt(synthss_);

		//// MATCHING TO ASSIGN CAPACTITY CONSTRAINTS

		//scene_.set_rt(synthss_);
	};

	return false;
}

bool CMeshSynthesizer::nnf_matching()
{
	// pre_compute
	compute_local_patches(synthss_, query_patches_, false);
	std::vector<CLocalPatch> best_nn_list;
	// do mrf_based_weight_assignment
	for (int i = 0; i < query_patches_.size(); i++)
	{
		CLocalPatch bnn;
		mrf_based_capacity_assignment(query_patches_[i], bnn);
		best_nn_list.push_back(bnn);
	}

	return false;
}

//void CMeshSynthesizer::direction_parallel_transport()
//{
//	std::cout << "to be constructed" << std::endl;
//}

void CMeshSynthesizer::compute_local_patches(Regular_triangulation* t, std::vector<CLocalPatch> &patch_list, bool is_example)
{
	std::cout << "compute local patches" << std::endl;

	patch_list.clear();
	for (auto vit = t->finite_vertices_begin(); vit != t->finite_vertices_end(); ++vit)
	{
		// skip boundary vertices
		bool is_boundary = false;
		std::vector<Vector_2> nn_vecs;
		auto vv_cir = t->incident_vertices(vit);
		auto vv_end = vv_cir;
		do
		{
			Vector_2 v = vv_cir->point().point() - vit->point().point();
			nn_vecs.push_back(v);
			if (t->is_infinite(vv_cir))
			{
				is_boundary = true;
				break;
			}
		} while (++vv_cir != vv_end);
		if (is_boundary)
			continue;

		CLocalPatch tmp_lp;
		tmp_lp.center_vh = vit;
		tmp_lp.is_example = is_example;
		tmp_lp.nn_vecs = nn_vecs;
		patch_list.push_back(tmp_lp);
	}
}

void CMeshSynthesizer::get_candidate_patch(CLocalPatch &bnn, int patch_id, int start_id)
{
	CLocalPatch elp = candidate_patches_[patch_id];
	elp.nn_shift(start_id);
	bnn = elp;
}

void CMeshSynthesizer::mrf_based_capacity_assignment(CLocalPatch &slp, CLocalPatch &bnn)
{
	CLocalPatch elp;
	
	int min_patch = -1;
	int min_start = -1;
	double min_error = std::numeric_limits<double>::max();

	for (int i = 0; i < candidate_patches_.size(); i++) 
	{
		elp = candidate_patches_[i];
		double errr = 0;
		int start_id = -1;
		errr = nnf_search(slp, elp, start_id);
		if (errr < min_error)
		{
			min_patch = i;
			min_start = start_id;
			min_error = errr;
		}
	}
	//
	get_candidate_patch(bnn, min_patch, min_start);

	std::cout << "assign weights... to be constructed..." << std::endl;
}

double CMeshSynthesizer::nnf_search(const CLocalPatch &slp, CLocalPatch &elp, int &start_id)
{
	int min_id = -1;
	double min_error = std::numeric_limits<double>::max();
	int end = elp.nn_vecs.size();
	for (int i = 0; i < end; i++)
	{
		// starting at i
		int iter_id = i;
		auto tmp_lp = elp;
		tmp_lp.nn_shift(iter_id);
		
		// skip 
		if (slp.nn_vecs.size() > tmp_lp.nn_vecs.size())
			continue;

		// else
		assert(slp.nn_vecs.size() <= tmp_lp.nn_vecs.size());
		double tmp_error = compute_matching_error(slp.nn_vecs, tmp_lp.nn_vecs);

		// record
		if (min_error < tmp_error)
		{
			min_id = iter_id;
			min_error = tmp_error;
		}
	}
	start_id = min_id;
	return min_error;
}

double CMeshSynthesizer::compute_matching_error(std::vector<Vector_2> s_nn_pts, std::vector<Vector_2> e_nn_pts)
{
	if (false)
	{
		std::cout << "use SVD error" << std::endl;
		return compute_svd_matching_error(s_nn_pts, e_nn_pts);
	}
	else
	{
		std::cout << "use matching error" << std::endl;
		return compute_index_matching_error(s_nn_pts, e_nn_pts);
	}
}

double CMeshSynthesizer::compute_svd_matching_error(std::vector<Vector_2> s_nn_pts, std::vector<Vector_2> e_nn_pts)
{
	assert(s_nn_pts.size() <= e_nn_pts.size());

	Eigen::MatrixXd snn_mat, enn_mat;
	snn_mat.resize(s_nn_pts.size(), 2);
	snn_mat.setZero();
	enn_mat.resize(e_nn_pts.size(), 2);
	enn_mat.setZero();

	for (int i = 0; i < s_nn_pts.size(); i++)
	{
		snn_mat(i, 0) = s_nn_pts[i].x();
		snn_mat(i, 1) = s_nn_pts[i].y();
	}
	for (int i = 0; i < e_nn_pts.size(); i++)
	{
		enn_mat(i, 0) = e_nn_pts[i].x();
		enn_mat(i, 1) = e_nn_pts[i].y();
	}

	//
	auto Res = CGeoCalculator::compute_rotation_from_right_to_left(snn_mat, enn_mat);
	auto rot_enn_mat = Res*enn_mat;
	double d = .0;
	for (int i = 0; i < s_nn_pts.size(); i++)
	{
		d += (s_nn_pts[i] - e_nn_pts[i]).squared_length();
	}
	d /= double(s_nn_pts.size());
	return d;
}

double CMeshSynthesizer::compute_index_matching_error(std::vector<Vector_2> s_nn_pts, std::vector<Vector_2> e_nn_pts)
{
	std::cout << "to be constructed" << std::endl;
	assert(s_nn_pts.size() <= e_nn_pts.size());
	double d = .0;
	for (int i = 0; i < s_nn_pts.size(); i++)
	{
		d += (s_nn_pts[i] - e_nn_pts[i]).squared_length();
	}
	d /= double(s_nn_pts.size());
	return d;
}