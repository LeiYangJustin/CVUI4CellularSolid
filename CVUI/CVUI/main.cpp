#include  <string>
#include "cv_window.h"
#include "../DataColle/voronoi_diagram.h"
#include "../DataColle/skeleton_extractor.h"
#include "../DataColle/reconstructor.h"
#include "../DataColle/mesh_optimizer.h"

int main(int argc, char** argv)
{
	// pre-defined parameters
	double Threshold = 1.0;
	double error = 100.0;

	/*----------------------------------------------*/
	// Step 1: extraction of skeletons
	// read some gray-scale image from the given path
	std::string filename;
	cv::Mat src_img = cv::imread(filename);
	// extract skeletons
	CSkeletonExtractor skeletonExtractor(src_img);
	std::vector<iPoint2> solid_skeleton = skeletonExtractor.getSolidSkeleton();
	/*----------------------------------------------*/

	/*----------------------------------------------*/
	//Step 2:
	//initialize a set of X from set Q
	std::vector<iPoint2> X;
	skeletonExtractor.getVoidSkeletonSamples(X);
	/*----------------------------------------------*/
	
	/*----------------------------------------------*/
	//Step 3: Two-step Optimization Loop
	do 
	{
		//Step 3.0: D <-- VoronoiDecomposer(X, Domain)
		//Regular_triangulation * rt = new Regular_triangulation;
		//rt->insert(X.begin(), X.end());

		CVoronoiDiagram *vd;
		vd->SetSites(X);

		//Step 3.1: D <-- Reconstructor(P, D)
		CReconstructor reconstructor;
		reconstructor.Update(vd);

		//Step 3.2: X <-- MeshOptimizer(D, X)
		CMeshOptimizer meshOptimizer;
		std::vector<iPoint2> Xnew;
		meshOptimizer.Update(vd, Xnew);
		
		//Step 3.3: compute error
		for (int i = 0; i < X.size(); i++)
		{
			Point<int> tmpP = Xnew[i] - X[i];
			error += tmpP.norm();
		}
		error /= double(X.size());

		// Set Xnew to X
		X = Xnew;

	} while (error > Threshold);
	




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

	//return 0;

}