#pragma once

#include <vector>
#include <iostream>

// include Eigen
#include "opt_prereq.h"
#include "my_utils.h"
#include "../LpCVT/lpcvt.h"

class OPT_CLASS CMyOptimizer
{
	//
	typedef Geex::vec3 vec3;

public:
	CMyOptimizer();
	~CMyOptimizer();

	//////////////
	// set info //
	//////////////

	void set_model(CLpCVT * cvt_obj) {
		cvt_object_ = cvt_obj;
	};

	void set_verbose_lvl(unsigned int verbose_lvl)
	{
		verbose_level_ = verbose_lvl;
	}

	//////////
	// loop //
	//////////
	
	void loop(unsigned int max_iter, double end_threshold)
	{
		unsigned int cnt_iter = 0;
		do
		{
			double norm = iterate(cnt_iter);
			std::cout << "current norm " << norm << " and end_threshold " << end_threshold << std::endl;
			std::cout << std::endl;
			if (norm < end_threshold)
			{
				std::cout << "reach local optimum" << std::endl;
				break;
			}
		} while (++cnt_iter < max_iter);
	};

private:

	/////////////
	// iterate //
	/////////////

	double iterate(int iter);

	////////////////
	// evaluation //
	////////////////
	
	// evaluation objective function
	double compute_objective();

	// evaluaton gradient of the objective function
	double compute_gradient(std::vector<vec3> &gradients);

	/////////////////
	// line search //
	/////////////////

	double do_line_search(const std::vector<vec3> &x0, const std::vector<vec3> &v0);

	////////////
	// update //
	////////////

	// update cvt_object_ inside
	void update_positions(
		std::vector<vec3> &x,
		const std::vector<vec3> &v, 
		const double step);

	void save_cvt_results() {
		is_export_model_ = true;
	}

private:
	CLpCVT * cvt_object_;
	unsigned int verbose_level_ = 1;
	bool is_export_model_ = false;
};

