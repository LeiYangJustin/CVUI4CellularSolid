#include <vector>
#include <iostream>

// INCLUDE OPENCV DIR
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include <opencv2/opencv.hpp>

//
#include <CGAL/bounding_box.h>

#define TEST_MAIN
#define THRESHOLD 1.0
#define MAXITER 3

//
#include "../DataColle/img_data.h"
#include "../AlgColle/skeleton_extractor.h"

//#include "../DataColle/voronoi_diagram.h"
//#include "../AlgColle/reconstructor.h"
#include "../AlgColle/mesh_optimizer.h"

//#include "../AlgColle/mesh_solver.h"
//#include "cv_window.h"
#include "reader.h"
#include "cv_viewer.h"

#ifdef TEST_MAIN

int main()
{
	// read some gray-scale image from the given path
	std::string filename = "D:\\MyProjects\\CVUI4CellularSolid\\CVUI\\img_data\\example2.png";
	cv::Mat src_img = cv::imread(filename);
	if (src_img.data == NULL)
		return EXIT_FAILURE;
	
#ifdef _DEBUG
	int height = 201;
	std::cout << "We are in the debugging mode" << std::endl;
#else 
	int height = 401;
#endif

	height = 201;
	int width = double(height) / double(src_img.rows)*double(src_img.cols);
	cv::resize(src_img, src_img, cv::Size(width, height));


	
	CImgData* img_data = new CImgData(src_img);

	std::cout << "\n-compute distance field from image" << std::endl;
	int rows, cols;
	std::vector<double> solid_field, void_field;
	img_data->get_two_distance_transform_fields(rows, cols, solid_field, void_field);
	//img_data->get_two_distance_transform_fields(rows, cols, void_field, solid_field);


	std::cout << "\n-optimize the triangulation" << std::endl;
	CMeshOptimizer mesh_optimizer;
	mesh_optimizer.background_setup(rows, cols, solid_field, void_field);
 	//mesh_optimizer.run_optimize_single_pt_as_test();
	mesh_optimizer.run_optimize();

	std::cout << "\n-output data" << std::endl;
}

//
//int main()
//{
//	// read some gray-scale image from the given path
//	std::string filename = "D:\\MyProjects\\CVUI4CellularSolid\\CVUI\\img_data\\example7.png";
//	cv::Mat src_img = cv::imread(filename);
//	if (src_img.data == NULL)
//		return EXIT_FAILURE;
//
//#ifdef _DEBUG
//	int height = 201;
//	std::cout << "We are in the debugging mode" << std::endl;
//#else 
//	int height = 401;
//#endif
//	int width = double(height) / double(src_img.rows)*double(src_img.cols);
//	cv::resize(src_img, src_img, cv::Size(width, height));
//	CImgData* img_data = new CImgData(src_img);
//
//	// Extraction of skeletons
//	CSkeletonExtractor skeletonExtractor(img_data);
//	skeletonExtractor.getSolidSkeleton();
//	skeletonExtractor.getVoidSkeleton();
//	skeletonExtractor.getImgData(img_data);
//
//	// Initialization of the voronoi sites
//	std::vector<cv::Point> samples;
//	skeletonExtractor.getVoidSkeletonSamples(samples);
//	//std::vector<Point_2> voro_sites;
//	std::vector<WPoint> voro_sites;
//	for (int i = 0; i < samples.size(); i++)
//	{
//		//voro_sites.push_back(Point_2(samples[i].x, samples[i].y));
//		WPoint wp(Point_2(samples[i].x, samples[i].y), 1.0);
//		voro_sites.push_back(wp);
//	}
//	std::set<Point_2> solid_skeletal_pts;
//	for (int i = 0; i < img_data->GetSolidSkeleton().size(); i++)
//	{
//		Point_2 p(img_data->GetSolidSkeleton()[i].x, img_data->GetSolidSkeleton()[i].y);
//		solid_skeletal_pts.insert(p);
//	}
//	std::set<Point_2> void_skeletal_pts;
//	for (int i = 0; i < img_data->GetVoidSkeleton().size(); i++)
//	{
//		Point_2 q(img_data->GetVoidSkeleton()[i].x, img_data->GetVoidSkeleton()[i].y);
//		void_skeletal_pts.insert(q);
//	}
//	Regular_triangulation *rt = new Regular_triangulation;
//	rt->insert(voro_sites.begin(), voro_sites.end());
//
//	CMeshSolver mesh_solver;
//	double tol = 0.5;
//	int max_iter = 200;
//	int every_x_step_to_update = 1;
//	mesh_solver.optimize(solid_skeletal_pts, void_skeletal_pts, rt, tol, max_iter, every_x_step_to_update);
//
//	CReader::write_triangles_from_rt("data.tria", rt);
//	CReader::write_faces_nearest("data.face_near", rt);
//	CReader::write_vertices_nearest("data.vert_near", rt);
//	CReader::write_pts("data.solid", solid_skeletal_pts);
//	CReader::write_pts("data.void", void_skeletal_pts);
//
//	delete rt;
//
//	return EXIT_SUCCESS;
//}

#else

int main()
{
	// read some gray-scale image from the given path
	std::string filename = "D:\\MyProjects\\CVUI4CellularSolid\\CVUI\\img_data\\test_example.png";
	cv::Mat src_img = cv::imread(filename);
	if (src_img.data == NULL)
		return EXIT_FAILURE;

#ifdef _DEBUG
	int height = 201;
	std::cout << "We are in the debugging mode" << std::endl;
#else 
	int height = 401;
#endif
	int width = double(height) / double(src_img.rows)*double(src_img.cols);
	cv::resize(src_img, src_img, cv::Size(width, height));
	CImgData* img_data = new CImgData(src_img);

	// Extraction of skeletons
	CSkeletonExtractor skeletonExtractor(img_data);
	skeletonExtractor.getSolidSkeleton();
	skeletonExtractor.getVoidSkeleton();
	skeletonExtractor.getImgData(img_data);

	// Initialization of the voronoi sites
	std::vector<cv::Point> samples;
	skeletonExtractor.getVoidSkeletonSamples(samples);
	std::vector<WPoint> voro_seeds;
	for (int i = 0; i < samples.size(); i++)
	{
		WPoint wp(samples[i].x, samples[i].y);
		voro_seeds.push_back(wp);
	}
	std::vector<Point_2> solid_skeletal_pts;
	for (int i = 0; i < img_data->GetSolidSkeleton().size(); i++)
	{
		Point_2 p(img_data->GetSolidSkeleton()[i].x, img_data->GetSolidSkeleton()[i].y);
		solid_skeletal_pts.push_back(p);
	}
	std::vector<Point_2> void_skeletal_pts;
	for (int i = 0; i < img_data->GetVoidSkeleton().size(); i++)
	{
		Point_2 q(img_data->GetVoidSkeleton()[i].x, img_data->GetVoidSkeleton()[i].y);
		void_skeletal_pts.push_back(q);
	}
	// constraint label for boundary handling
	std::vector<bool> is_constrained_list(voro_seeds.size(), false);
	
	//Step 3: Two-step Optimization Loop	
	std::cout << "number of pts to be reconstructed " << solid_skeletal_pts.size() << std::endl;
	CReconstructor reconstructor;
	reconstructor.SetReconstructionPts(solid_skeletal_pts);

	std::cout << "number of pts as constraint " << void_skeletal_pts.size() << std::endl;
	CMeshOptimizer meshOptimizer;
	meshOptimizer.SetConstraintPts(void_skeletal_pts);

	std::vector<std::vector<Segment_2>> history_tria_edges;
	std::vector<std::vector<Segment_2>> history_voro_edges;
	int cntIter = 0;

	CVoronoiDiagram* pVD = new CVoronoiDiagram;
	pVD->SetBoundingBox(width, height);
	//
	do 
	{
		//Step 3.0: D <-- VoronoiDecomposer(X, Domain)
		std::cout << "Build voronoi diagram with " << voro_seeds.size() << " of voronoi sites" << std::endl;
		pVD->SetSeedSites(voro_seeds, is_constrained_list);

		// data for plot
		std::vector<Segment_2> tria_edges, voro_edges;
		pVD->GetCroppedTriangulatedSegments(tria_edges);
		std::cout << tria_edges.size() << std::endl;
		pVD->GetCroppedVoronoiSegments(voro_edges);
		std::cout << voro_edges.size() << std::endl;
		history_tria_edges.push_back(tria_edges);
		history_voro_edges.push_back(voro_edges);

		//Step 3.1: PD <-- Reconstructor(S1, D)
		std::cout << "Reconstruction" << std::endl;
		reconstructor.Update(pVD);
		std::cout << "Reconstruction error: " << reconstructor.GetReconstructionError() << std::endl;

		//// data for plot
		//std::vector<Point_2> rp_list;
		//std::vector<int> rp_label_list;
		//reconstructor.GetReconstructionPts(rp_list, rp_label_list);

		//Step 3.2: Xnew <-- MeshOptimizer(D, X)
		std::cout << "Mesh Optimization" << std::endl;
		std::vector<WPoint> old = voro_seeds;
		meshOptimizer.Update2(pVD, voro_seeds);

		double error = 0.0;
		for (int i = 0; i < voro_seeds.size(); i++)
		{
			Vector_2 tmpP = voro_seeds[i].point() - old[i].point();
			double tmpW = voro_seeds[i].weight() - old[i].weight();
			error += (tmpP.squared_length()) /*+ tmpW*tmpW*/;
			//std::cout << tmpW << std::endl;
		}
		error /= double(voro_seeds.size());
		std::cout << "Loop = " << cntIter << "; averaged error = " << error << std::endl;
		std::cout << std::endl;

	} while (/*error > Threshold &&*/ ++cntIter < MAXITER);
	std::cout << "We are done here" << std::endl;

	// PLOT
	CViewer viewer(*img_data);
	for (int i = 0; i < cntIter; i++)
	{
		viewer.clearData();

		//std::cout << "tria edges: " << std::endl;
		for (auto edge : history_tria_edges[i])
		{
			viewer.addSegment(edge, 2);
			//std::cout << edge << std::endl;
		}
		//std::cout << std::endl;
		//std::cout << "voronoi edges: " << std::endl;
		for (auto edge : history_voro_edges[i])
		{
			viewer.addSegment(edge, 7);
			//std::cout << edge << std::endl;
		}
		//std::cout << std::endl;
		std::string winName = "Test";
		winName.append(std::to_string(i));
		viewer.draw(winName);
	}

	// end
	//delete pVD;
	delete img_data;

	return EXIT_SUCCESS;
}

#endif // TEST_MAIN