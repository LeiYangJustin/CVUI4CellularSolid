#include <vector>
#include <iostream>

// INCLUDE OPENCV DIR
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include <opencv2/opencv.hpp>

//
#include "../DataColle/img_data.h"
#include "../DataColle/voronoi_diagram.h"
#include "../AlgColle/skeleton_extractor.h"
#include "../AlgColle/reconstructor.h"
#include "../AlgColle/mesh_optimizer.h"
#include "cv_window.h"

// pre-defined parameters
double Threshold = 1.0;
double error = 100.0;
int cntIter = 0;
int MaxIter = 1;

int main()
{
	// read some gray-scale image from the given path
	std::string filename = "D:\\MyProjects\\CVUI4CellularSolid\\CVUI\\img_data\\test_example.png";
	cv::Mat src_img = cv::imread(filename);
	if (src_img.data == NULL)
		return EXIT_FAILURE;

#ifdef _DEBUG
	int height = 401;
	std::cout << "We are in the debugging mode" << std::endl;
#else 
	int height = 401;
#endif
	// resize the image
	int width = double(height) / double(src_img.rows)*double(src_img.cols);
	cv::resize(src_img, src_img, cv::Size(width, height));

	// set image data
	CImgData* img_data = new CImgData(src_img);

	/*----------------------------------------------*/
	// Step 1: extraction of skeletons
	CSkeletonExtractor skeletonExtractor(img_data);
	skeletonExtractor.getSolidSkeleton();
	skeletonExtractor.getVoidSkeleton();
	skeletonExtractor.getImgData(img_data);
	/*----------------------------------------------*/

	/*----------------------------------------------*/
	//Step 2:
	//initialize a set of X from set Q
	std::vector<cv::Point> samples;
	skeletonExtractor.getVoidSkeletonSamples(samples);
	/*----------------------------------------------*/


	// Data conversion from OPENCV to CGAL
	std::vector<WPoint> X;
	for (int i = 0; i < samples.size(); i++)
	{
		WPoint wp(samples[i].x, samples[i].y);
		X.push_back(wp);
	}
	std::vector<Point_2> P;
	for (int i = 0; i < img_data->GetSolidSkeleton().size(); i++)
	{
		Point_2 p(img_data->GetSolidSkeleton()[i].x, img_data->GetSolidSkeleton()[i].y);
		P.push_back(p);
	}
	std::vector<Point_2> Q;
	for (int i = 0; i < img_data->GetVoidSkeleton().size(); i++)
	{
		Point_2 q(img_data->GetVoidSkeleton()[i].x, img_data->GetVoidSkeleton()[i].y);
		Q.push_back(q);
	}

	/*----------------------------------------------*/
	//Step 3: Two-step Optimization Loop
	// declaration of the components
	CVoronoiDiagram* pVD = new CVoronoiDiagram;
	pVD->UpdateTriangulation(X);

	CReconstructor reconstructor;
	reconstructor.SetBBox(width, height);

	std::cout << "number of pts to be reconstructed " << P.size() << std::endl;
	reconstructor.SetReconstructionPts(P);

	CMeshOptimizer meshOptimizer;
	meshOptimizer.SetConstraintPts(Q);

	do {
		std::cout << "#VoroSites: " << X.size() << std::endl;
		
		std::cout << "UpdateTriangulation" << std::endl;
		//Step 3.0: D <-- VoronoiDecomposer(X, Domain)
		pVD->UpdateTriangulation(X);
		
		std::cout << "Reconstruction" << std::endl;
		//Step 3.1: D <-- Reconstructor(P, D)
		reconstructor.Update(pVD);
		// get reconstruction accuracy
		std::cout << "Reconstruction error: " << reconstructor.GetReconstructionError() << std::endl;

		std::cout << "Mesh Optimization" << std::endl;
		//Step 3.2: Xnew <-- MeshOptimizer(D, X)
		std::vector<WPoint> Xold = X;
		meshOptimizer.Update(pVD, X);


		for (int i = 0; i < X.size(); i++)
		{
			Vector_2 tmpP = Xold[i].point() - X[i].point();
			double tmpW = Xold[i].weight() - X[i].weight();
			error += (tmpP.squared_length() + tmpW*tmpW);
		}
		error /= double(X.size());
		std::cout << cntIter << std::endl;
	} while (error > Threshold || ++cntIter < MaxIter);
	// Repeat until criteria met
	std::cout << "We are done here" << std::endl;



	// Test
	std::vector<Point_2> rp_list;
	std::vector<int> rp_label_list;
	reconstructor.GetReconstructionPts(rp_list, rp_label_list);

	// PLOT
	CVoronoiDrawer* vd_drawer = new CVoronoiDrawer;
	vd_drawer->SetImgData(img_data);
	vd_drawer->SetVD(pVD);
	CCVWindow cvWin;
	cvWin.SetVoronoiDrawer(vd_drawer);
	cvWin.ShowReconstructionPointClusters(rp_list, rp_label_list, 0);
	cvWin.ShowUpdatedVoronoiDiagram(0);


	return EXIT_SUCCESS;
}