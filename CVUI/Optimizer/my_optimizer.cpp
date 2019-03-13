#include "my_optimizer.h"
#include "my_line_search.h"

CMyOptimizer::CMyOptimizer()
{
}

CMyOptimizer::~CMyOptimizer()
{
}

double CMyOptimizer::iterate(int iter)
{
	std::cout << iter << " iteration..." << std::endl;

	save_cvt_results();

	// collect positions
	std::cout << "--collect sites..." << std::endl;
	std::vector<vec3> points;
	cvt_object_->get_sites(points);

	// compute gradient
	std::cout << "--compute gradient..." << std::endl;
	std::vector<vec3> gradients;
	double gnorm = compute_gradient(gradients);

	std::vector<vec3> negative_gradients;
	for (unsigned int i = 0; i < gradients.size(); ++i)
		negative_gradients.push_back(-1 * gradients[i]);

	// line search
	std::cout << "--do line search..." << std::endl;
	double step = do_line_search(points, negative_gradients);
	std::cout << "Line_search done: step = " << step << std::endl;

	// update
	std::cout << "--update sites..." << std::endl;
	update_positions(points, negative_gradients, step);

	cvt_object_->update_decomposition(points, is_export_model_);

	// do we need to recompute the gradient
	return gnorm;
};

double CMyOptimizer::compute_objective()
{
	double fval = cvt_object_->cvt_energy();
	return fval;
}

double CMyOptimizer::compute_gradient(std::vector<vec3> &gradients)
{
	gradients.clear();
	cvt_object_->cvt_gradient(gradients);
	double gnorm = Utils::normalize_gradients(gradients);
	return gnorm;
}

double CMyOptimizer::do_line_search(const std::vector<vec3> &x0,  const std::vector<vec3> &v0)
{
	const int max_ls_iter = 10;
	double max_alpha = 1.0;
	MyLineSearch line_searcher(cvt_object_, max_ls_iter, max_alpha);
	double step = line_searcher.search_with_backtracking(x0, v0);
	return step;
}

void CMyOptimizer::update_positions(std::vector<vec3> &x, const std::vector<vec3> &v, const double step)
{
	assert(x.size() == v.size());

	//
	std::vector<vec3> xx(x.size());
	for (unsigned int i = 0; i < x.size(); ++i)
		xx[i] = x[i] + step*v[i];

	x = xx;
}
