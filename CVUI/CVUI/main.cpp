#include  <string>
#include "cv_window.h"
#include "../DataColle/voronoi_diagram.h"
#include "../AlgColle/skeleton_extractor.h"
#include "../AlgColle/reconstructor.h"
#include "../AlgColle/mesh_optimizer.h"

int main(int argc, char** argv)
{
	std::cout << "hello world" << std::endl;

	// pre-defined parameters
	double Threshold = 1.0;
	double error = 100.0;
	int cntIter = 0;
	int MaxIter = 1;

	/*----------------------------------------------*/
	// Step 1: extraction of skeletons
	// read some gray-scale image from the given path
	std::string filename = "D:\\MyProjects\\CVUI4CellularSolid\\CVUI\\img_data\\test_example.png";
	cv::Mat src_img = cv::imread(filename);
	if (src_img.data == NULL)
		return 0;

	CImgData img_data;
	img_data.SetSrcImg(src_img);

	// extract skeletons
	CSkeletonExtractor skeletonExtractor(img_data);
	std::vector<cv::Point> solid_skeleton = skeletonExtractor.getSolidSkeleton();
	std::vector<cv::Point> void_skeleton = skeletonExtractor.getVoidSkeleton();
	/*----------------------------------------------*/

	/*----------------------------------------------*/
	//Step 2:
	//initialize a set of X from set Q
	std::vector<cv::Point> samples;
	skeletonExtractor.getVoidSkeletonSamples(samples);
	/*----------------------------------------------*/
	
	/*----------------------------------------------*/
	//Step 3: Two-step Optimization Loop
	// Data conversion from OPENCV to CGAL
	std::vector<WPoint> X, Xnew;
	for (int i = 0; i < samples.size(); i++)
	{
		WPoint wp(samples[i].x, samples[i].y);
		Xnew.push_back(wp);
	}
	std::vector<Point2> P;
	for (int i = 0; i < solid_skeleton.size(); i++)
	{
		Point2 p(solid_skeleton[i].x, solid_skeleton[i].y);
		P.push_back(p);
	}
	// declaration of the components
	CReconstructor reconstructor;
	reconstructor.SetReconstructionPts(P);
	CMeshOptimizer meshOptimizer;
	CVoronoiDiagram *vd = new CVoronoiDiagram;
	do 
	{
		// 
		X = Xnew;
		std::cout << "#VoroSites: " << X.size() << std::endl;
		//Step 3.0: D <-- VoronoiDecomposer(X, Domain)
		vd->UpdateTriangulation(X);
		std::cout << "UpdateTriangulation" << std::endl;
		//Step 3.1: D <-- Reconstructor(P, D)
		reconstructor.Update(vd);
		std::cout << "Reconstruction" << std::endl;

		//Step 3.2: Xnew <-- MeshOptimizer(D, X)
		std::vector<WPoint> Xnew;
		meshOptimizer.Update(vd, Xnew);
		
		//Step 3.3: compute error
		// get reconstruction accuracy
		reconstructor.GetAccuracy();

		for (int i = 0; i < X.size(); i++)
		{
			Vector2 tmpP = Xnew[i].point() - X[i].point();
			double tmpW = Xnew[i].weight() - X[i].weight();
			error += (tmpP.squared_length() + tmpW*tmpW);
		}
		error /= double(X.size());	
		std::cout << cntIter << std::endl;
	} while (error > Threshold || ++cntIter < MaxIter);
	// Repeat until criteria met

	std::cout << "We are done here" << std::endl;

	//CCVWindow cvWindow(argv[0]);

	//// user interface
	//while (true)
	//{
	//	// add points
	//	int x = 0, y = 0;
	//	int flag;
	//	// get x, y and flag from the CV window
	//	flag = cvWindow.getMouseClick(x, y);
	//	cvWindow.showUpdatedVoronoiDiagram();
	//}



	delete vd;
	//return 0;

}