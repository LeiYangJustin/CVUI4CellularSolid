#include <iostream>

#include "../Triangulation3/tria3_optimizer.h"
#include "../Optimizer/my_optimizer.h"
#include "lpcvt_wrapper.h"

int main()
{
	std::cout << "3D MODELING" << std::endl;

	//// Filename
	//std::string mesh_filename = "D:\\MyProjects\\CVUI4CellularSolid\\CVUI\\MainPrj3D\\test_data\\three_holes.obj";
	//std::string pts_filename = "D:\\MyProjects\\CVUI4CellularSolid\\CVUI\\MainPrj3D\\test_data\\three_holes.pts";
	std::string mesh_filename = "D:\\MyProjects\\CVUI4CellularSolid\\CVUI\\MainPrj3D\\test_data\\cube.obj";

	// parameters
	unsigned int num_sites = 100;
	unsigned int max_iter = 100;
	unsigned int p_order = 2;
	double end_threshold = 0.001;
	MyMesh mesh;
	OpenMesh::IO::read_mesh(mesh, mesh_filename);

	// new wrapper for optimization
	CLPCVTwrapper* cvt_wrapper = new CLPCVTwrapper;
	cvt_wrapper->init(num_sites, mesh, p_order);
	cvt_wrapper->optimize(max_iter, end_threshold);

	// end
	system("pause");
}