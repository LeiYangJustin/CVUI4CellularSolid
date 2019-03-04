#ifndef C_EXTRACTION_H
#define C_EXTRACTION_H

#include "algprereq.h"

#include <vector>

//
#include <Eigen/Dense>
#include <Eigen/Sparse>

// INCLUDE OPENCV DIR
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include <opencv2/opencv.hpp>

#include "../DataColle/types.h"
#include "bg_scene.h"

class ALGCOLLE_CLASS CMeshExtractor
{
public:
	CMeshExtractor(RT *rt) { rt_ = rt; };
	~CMeshExtractor() { delete rt_; };

public:
	void setup_background(int rows, int cols, 
		std::vector<double> &sz, std::vector<double> &vz, std::vector<std::vector<double>> edgePts);
	void run_extraction();

	// examination
	bool write_input_fields_for_check(Eigen::MatrixXd &z, Eigen::MatrixXd &zgx, Eigen::MatrixXd &zgy);

private:
	// access distance transform data
	void compute_field_gradients(
		const Eigen::MatrixXd &zmap,
		Eigen::MatrixXd &zgradx,
		Eigen::MatrixXd &zgrady);

	void compute_principal_directions(
		const Eigen::MatrixXd &zmap,
		Eigen::MatrixXd &pc1x,
		Eigen::MatrixXd &pc1y,
		Eigen::MatrixXd &pc2x,
		Eigen::MatrixXd &pc2y);

private:
	RT * rt_;

	int rows_;
	int cols_;

	Eigen::MatrixXd z_; // distance transform 
	Eigen::MatrixXd zgx_; // distance transform gradient
	Eigen::MatrixXd zgy_; // distance transform gradient

	//std::vector<double> x_grid_;
	//std::vector<double> y_grid_;

	BgScene bg_scene_;

	std::vector<double> feval_history_;

};

#endif // !C_EXTRACTION_H
