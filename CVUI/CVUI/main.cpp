#include <vector>
#include <iostream>

// INCLUDE OPENCV DIR
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include <opencv2/opencv.hpp>
#include <CGAL/bounding_box.h>

//#define TEST_EXAMPLE
//#define TEST_SYNTHSS
#define THRESHOLD 1.0
#define MAXITER 3

//
#include "../DataColle/img_data.h"
//#include "../AlgColle/skeleton_extractor.h"
#include "../AlgColle/dGBOD12_test.h"

#include "reader.h"
#include "cv_viewer.h"



//int main()
//{
//	// read some gray-scale image from the given path
//	std::string filename = "D:\\MyProjects\\CVUI4CellularSolid\\CVUI\\img_data\\example2.png";
//	cv::Mat src_img = cv::imread(filename);
//	if (src_img.data == NULL)
//		return EXIT_FAILURE;
//	
//	//
//#ifdef _DEBUG
//	int height = 201;
//	std::cout << "We are in the debugging mode" << std::endl;
//#else 
//	int height = 401;
//#endif
//	height = 201;
//	int width = double(height) / double(src_img.rows)*double(src_img.cols);
//	cv::resize(src_img, src_img, cv::Size(width, height));
//	CImgData* img_data = new CImgData(src_img);
//
//	//
//	std::cout << "\n-compute distance field from image" << std::endl;
//	int rows, cols;
//	std::vector<double> solid_field, void_field;
//	img_data->get_two_distance_transform_fields(rows, cols, solid_field, void_field);
//
//	//
//	std::cout << "\n-optimize the triangulation" << std::endl;
//	Regular_triangulation * example_net = new Regular_triangulation;
//
//#ifdef TEST_EXAMPLE
//	CMeshOptimizer mesh_optimizer(example_net);
//	mesh_optimizer.background_setup(rows, cols, solid_field, void_field);
//	mesh_optimizer.run_optimize();
//
//	std::cout << "\n-output data" << std::endl;
//	std::cout << "#verts = " << example_net->number_of_vertices() << std::endl;
//#endif
//
//#ifdef TEST_SYNTHSS
//	// 
//	std::cout << "\n-do synthesis..." << std::endl;
//	double height_ss = 100;
//	double width_ss = 100;
//
//	std::vector<Point_2> pts;
//	CReader::read_data_array_file("D:\\MyProjects\\CVUI4CellularSolid\\CVUI\\img_data\\samples_40.txt", pts);
//	Regular_triangulation * synthss_net = new Regular_triangulation;
//	CMeshSynthesizer mesh_synthesizer(synthss_net);
//
//#endif
//
//	//std::cout << "finished" << std::endl;
//	//example_net->clear();
//	//synthss_net->clear();
//	//delete example_net;
//	//delete synthss_net;
//}


#include "../AlgColle/synthesis.h"

int main()
{
	//std::vector<Point_2> pts;
	//CReader::read_data_array_file("D:\\MyProjects\\CVUI4CellularSolid\\CVUI\\img_data\\samples_40.txt", pts);

	//CTest_dGBOD12 test;
	//test.run_test(200);

	RT * example = new RT;
	Domain in_domain;
	in_domain.init_rectangle(1, 1);
	in_domain.init_area();

	CMeshSynthesizer synthesizer;
	synthesizer.set_exemplar(example);
	synthesizer.set_domain(in_domain);
	synthesizer.synthesize();



}