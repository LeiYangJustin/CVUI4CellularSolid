#ifndef _BG_SCENE_H_
#define _BG_SCENE_H_

// STL
#include <vector>

// EIGEN
#include <Eigen/Dense>
#include <Eigen/Sparse>

// local
#include "../DataColle/types.h"
#include "../DataColle/util.h"
#include "bbnot/scene.h"
#include "bbnot/pw_line_search.h"

#include "algprereq.h"

//// build a kdtree for spatial query
#include "ann_tree_wrapper.h"

//typedef CLSWeights<Scene, FT> LSWeights;
class BgScene;
typedef C_BG_LSPositions<BgScene, Point_2, Vector_2> BG_LSPositions;

class BgScene : public Scene
{
private:
	double alpha_;
	double beta_;

	std::vector<double> alm_lambda_;
	double alm_miu_;

	double rows_;
	double cols_;

	std::map<Point_2, FT> z_map_;
	std::map<Point_2, FT> zgx_map_;
	std::map<Point_2, FT> zgy_map_;

	bool has_edgepts_kdtree_ = false;
	ANNTreeWrapper tree_search_edgepts_tool_;

	std::vector<double> x_interval_;
	std::vector<double> y_interval_;

	//Eigen::MatrixXd z_; // distance transform 
	//Eigen::MatrixXd zgx_; // distance transform gradient
	//Eigen::MatrixXd zgy_; // distance transform gradient
	//std::vector<double> x_grid_;
	//std::vector<double> y_grid_;

	enum SolverType { CVT_SOLVER, FIELD_SOLVER, COMBINED_CVT_FIELD_SOLVER, ALM_SOLVER };
	SolverType solver_type_;

public:
	BgScene() : Scene() {};

	~BgScene()
	{
		clear();
	}
	
	void construct_bg_info(double cols, double rows,
		Eigen::MatrixXd &z, Eigen::MatrixXd &zgx, Eigen::MatrixXd &zgy,
		std::vector<Point_2> edge_pts);
	double alpha() { return alpha_; };
	double beta() { return beta_; };
	SolverType get_solvertype() const { return solver_type_; };

	// interface to set a user-prescribed domain
	void set_domain()
	{
		m_domain.clear();
		m_domain.init_rectangle(cols_, rows_);
		m_domain.init_area();
		m_rt.set_domain(&m_domain);
	};

	unsigned bg_optimize_all_for_extraction(FT wstep, FT xstep,
		unsigned max_newton_iters,
		FT epsilon, unsigned max_iters,
		std::ostream& out);

	// optimize with specific constraint
	FT bg_optimize_positions_via_combined_gradient_ascent(FT timestep,
		bool update);
	FT bg_optimize_positions_via_lagrangian_multiplier();

	// basic components
	void init_alm_parameters();

	// gradient solvers	
	void bg_compute_position_gradient(std::vector<Vector_2>& gradient, const FT coef, SolverType gstype);
	void bg_compute_position_gradient_cvt_solver(std::vector<Vector_2>& gradient, const FT coef);
	void bg_compute_position_gradient_field_solver(std::vector<Vector_2>& gradient, const FT coef);
	void bg_compute_position_gradient_combined_cvt_field_solver(std::vector<Vector_2>& gradient, const FT coef);
	void bg_compute_position_gradient_alm_solver(std::vector<Vector_2>& gradient, const FT coef, const std::vector<double> lambda);

	// energy solvers
	FT bg_compute_energy(SolverType stype);
	FT bg_compute_energy_wcvt_solver();
	FT bg_compute_energy_field_solver();
	FT bg_compute_energy_combined_cvt_field_solver();
	FT bg_compute_energy_alm_solver(const std::vector<double> lambda);
	

	double dist_transform_field_S(const Point_2 &pos)
	{
		return get_interp_val_from_2dgrid(pos, z_map_, 1.0);
	};
	double dist_transform_field_V(const Point_2 &pos)
	{
		return get_interp_val_from_2dgrid(pos, z_map_, -1.0);
	};
	Eigen::Vector2d dist_transform_grad_S(const Point_2 &pos)
	{
		//std::cout << "gx" << std::endl;
		double gx = get_interp_val_from_2dgrid(pos, zgx_map_, 1.0);
		//std::cout << "gy" << std::endl;
		double gy = get_interp_val_from_2dgrid(pos, zgy_map_, 1.0);
		return Eigen::Vector2d(gx, gy);
	};
	Eigen::Vector2d dist_transform_grad_V(const Point_2 &pos)
	{
		double gx = get_interp_val_from_2dgrid(pos, zgx_map_, -1.0);
		double gy = get_interp_val_from_2dgrid(pos, zgy_map_, -1.0);
		return Eigen::Vector2d(gx, gy);
	};
	void transform_by_scaled_eigen_vectors(Vector_2 &gradient, 
		const std::vector<Point_2> datapts, const Point_2 pt_cntr = Point_2());
	void project_gradient_to_feasible_direction(Vector_2 &g, const Point_2 p);

	//double get_interp_val_from_kdtree(const Point_2 &pos, std::map<Point_2, FT> &map, const double coef = 1.0);
	double get_interp_val_from_2dgrid(const Point_2 &pos, std::map<Point_2, FT> &map, const double coef = 1.0);
	void get_interval(double x, double &xl, double &xu, const std::vector<double> &interval);
	// post-processing
	void remove_seeds_in_proximity(std::vector<Point_2> &new_sites, double merging_length = 5);
	void insert_seeds_at_edges(std::vector<Point_2> &new_sites, double merging_length = 5);
	void insert_seeds_at_faces();

	// check
	void print_m_vertices() {
		for (int i = 0; i < m_vertices.size(); i++)
			std::cout << m_vertices[i]->point().point() << std::endl;
	}
};

#endif // _BG_SCENE_H_