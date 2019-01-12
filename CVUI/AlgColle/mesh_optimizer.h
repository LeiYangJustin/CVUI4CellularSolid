#ifndef C_MESH_OPTIMIZER_H
#define C_MESH_OPTIMIZER_H

#include "algprereq.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include "../DataColle/customized_tds.h"

class ALGCOLLE_CLASS CMeshOptimizer
{
public:
	CMeshOptimizer() { rt_ = new Regular_triangulation; };
	~CMeshOptimizer() { delete rt_; };

public:
	void background_setup(int rows, int cols, std::vector<double> sz, std::vector<double> vz);
	void run_optimize();
	void run_optimize_single_pt_as_test();

	// examination
	bool write_input_fields_for_check();
	bool write_updated_triangulation(std::string fname);
	bool write_weighted_vertices(std::string fname);
	bool write_updated_voronoi(std::string fname);
	bool write_single_pt_trace(std::string fname, std::vector<Point_2> trace);
	bool write_energy_history(std::string fname, std::vector<double> energy_history);

private:
	// removing duplicated points in wpts, construct a T, and insert new points if needed
	void rebuild_triangulation(std::vector<Point_2> wpts, double merging_length = 10);

	// optimization
	double update_vertex_i(Regular_triangulation::Vertex_handle vh, double init_stepsize,
		const double alpha, const double beta, const double gamma, bool use_cvt);
	double update_vertex_i_alternating(Regular_triangulation::Vertex_handle vh, double init_stepsize,
		const double alpha, const double beta, const double gamma, bool use_cvt);

	// function evaluation
	double evaluate_F(const Eigen::Vector2d &pos, Regular_triangulation::Vertex_handle vh, 
		const double alpha, const double beta, const double gamma, bool use_cvt);
	double evaluate_F_ODT(const Eigen::Vector2d &pos, Regular_triangulation::Vertex_handle vh,
		const double alpha, const double beta, const double gamma);
	double evaluate_F_CVT(const Eigen::Vector2d &pos, Regular_triangulation::Vertex_handle vh);
	double evaluate_F_constraint(const Eigen::Vector2d &pos, Regular_triangulation::Vertex_handle vh,
		const double alpha, const double beta);

	// gradient computation
	Eigen::Vector2d evaluate_G(const Eigen::Vector2d &pos, Regular_triangulation::Vertex_handle vh,
		const double alpha, const double beta, const double gamma, bool use_cvt);
	// using ODT energy + constraint
	Eigen::Vector2d evaluate_G_ODT(const Eigen::Vector2d &pos, Regular_triangulation::Vertex_handle vh, 
		const double alpha, const double beta, const double gamma);
	// compute only gradient of the CVT energy 
	Eigen::Vector2d evaluate_G_CVT(const Eigen::Vector2d &pos, Regular_triangulation::Vertex_handle vh);
	Eigen::Vector2d evaluate_G_constraint(const Eigen::Vector2d &pos, Regular_triangulation::Vertex_handle vh,
		const double alpha, const double beta);

	double find_feasible_step(const Eigen::Vector2d &pos, Regular_triangulation::Vertex_handle vh, const Eigen::Vector2d &dir);
	double find_steplength_with_wolfe_cond(double s,
		const Eigen::Vector2d &xi,
		const Eigen::Vector2d &gi,
		const Regular_triangulation::Vertex_handle vh,
		const double alpha, const double beta, const double gamma, bool use_cvt);

	//
	void insert_seeds_at_edges();
	void insert_seeds_at_faces();
	
	void boundary_sampling(Point_2 psrc, Point_2 ptgt, std::vector<Point_2> &samples);
	
	//// to be done
	//// amplitude prescribes the amplitude of the oscillation
	//void perturb_vertices(double amplitude);
	//// N prescribes the number of points being removed and added
	//void random_insert_and_remove(int N);

	// access distance transform data
	void CMeshOptimizer::compute_field_gradients(
		const Eigen::MatrixXd &zmap,
		Eigen::MatrixXd &zgradx,
		Eigen::MatrixXd &zgrady);

	double dist_transform_field_S(const Eigen::Vector2d &pos)
	{
		return get_interp_val_from_2dgrid(pos, z_);
	};
	double dist_transform_field_V(const Eigen::Vector2d &pos)
	{
		return get_interp_val_from_2dgrid(pos, -z_);
	};
	Eigen::Vector2d dist_transform_grad_S(const Eigen::Vector2d &pos)
	{
		double gx = get_interp_val_from_2dgrid(pos, zgx_);
		double gy = get_interp_val_from_2dgrid(pos, zgy_);
		return Eigen::Vector2d(gx, gy);
	};
	Eigen::Vector2d dist_transform_grad_V(const Eigen::Vector2d &pos)
	{
		double gx = get_interp_val_from_2dgrid(pos, -zgx_);
		double gy = get_interp_val_from_2dgrid(pos, -zgy_);
		return Eigen::Vector2d(gx, gy);
	};
	
	double get_interp_val_from_2dgrid(const Eigen::Vector2d &pos, const Eigen::MatrixXd &map);

private:
	Regular_triangulation * rt_;

	int rows_;
	int cols_;

	Eigen::MatrixXd z_; // distance transform 
	Eigen::MatrixXd zgx_; // distance transform gradient
	Eigen::MatrixXd zgy_; // distance transform gradient

	//Eigen::MatrixXd void_xy_z_; // distance transform 
	//Eigen::MatrixXd void_xy_gx_; // distance transform gradient
	//Eigen::MatrixXd void_xy_gy_; // distance transform gradient

	//Eigen::MatrixXd solid_xy_z_; // distance transform 
	//Eigen::MatrixXd solid_xy_gx_; // distance transform gradient
	//Eigen::MatrixXd solid_xy_gy_; // distance transform gradient

	std::vector<double> x_grid_;
	std::vector<double> y_grid_;

	std::vector<double> feval_history_;
};

#endif // !C_MESH_OPTIMIZER_H
