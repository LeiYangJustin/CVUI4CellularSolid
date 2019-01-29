#ifndef C_EXTRACTION_H
#define C_EXTRACTION_H

#include "algprereq.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>

#include "../DataColle/types.h"
#include "bbnot/scene.h"

class ALGCOLLE_CLASS CExtraction
{
public:
	CExtraction(RT *rt) { rt_ = rt; };
	~CExtraction() { delete rt_; };

public:
	void setup_background(int rows, int cols, std::vector<double> sz, std::vector<double> vz);
	void run_optimize();

private:
	//
	void remove_seeds_in_proximity(double merging_length = 1);
	void insert_seeds_at_edges();
	void insert_seeds_at_faces();

	// access distance transform data
	void compute_field_gradients(
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
	RT * rt_;

	int rows_;
	int cols_;

	Eigen::MatrixXd z_; // distance transform 
	Eigen::MatrixXd zgx_; // distance transform gradient
	Eigen::MatrixXd zgy_; // distance transform gradient

	std::vector<double> x_grid_;
	std::vector<double> y_grid_;

	Scene scene_;

	std::vector<double> feval_history_;
};

#endif // !C_EXTRACTION_H
