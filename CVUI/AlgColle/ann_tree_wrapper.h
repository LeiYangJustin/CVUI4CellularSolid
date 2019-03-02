#ifndef C_ANN_TREE_WRAPPER_H
#define C_ANN_TREE_WRAPPER_H

#include "algprereq.h"
#include "ANN.h"
#include <vector>
#include "../DataColle/types.h"

class ALGCOLLE_CLASS ANNTreeWrapper
{
public:
	ANNTreeWrapper(std::vector<Point_2> &dataPts) : eps(0) { build_tree(dataPts); };
	~ANNTreeWrapper() { delete kdTree; };
private:
	void build_tree(std::vector<Point_2> &dataPts);
	void get_k_nn_from_tree(Point_2 &in_queryPt, std::vector<Point_2> &NNPts, int k = 1);

	// Data
	int					nPts;					// actual number of data points
	ANNpointArray		dataPts;				// data points

	// Query
	double				eps;					// error bound

	// Tree
	ANNkd_tree*			kdTree;					// search structure
};

void ANNTreeWrapper::build_tree(std::vector<Point_2> &in_dataPts)
{
	int nPts = in_dataPts.size();
	int dim = 2;
	dataPts = annAllocPts(nPts, dim);			// allocate data points

	for (int i = 0; i < in_dataPts.size(); ++i)
	{
		dataPts[i][0] = in_dataPts[i].x();
		dataPts[i][1] = in_dataPts[i].y();
	}

	kdTree = new ANNkd_tree(		// build search structure
		dataPts,					// the data points
		nPts,						// number of points
		dim);						// dimension of space
}

void ANNTreeWrapper::get_k_nn_from_tree(Point_2 &in_queryPt, std::vector<Point_2> &NNPts, int k = 1)
{
	auto queryPt = annAllocPt(2);					// allocate query point
	auto nnIdx = new ANNidx[k];						// allocate near neigh indices
	auto dists = new ANNdist[k];					// allocate near neighbor dists

	queryPt[0] = in_queryPt.x();
	queryPt[1] = in_queryPt.y();

	kdTree->annkSearch(						// search
			queryPt,						// query point
			k,								// number of near neighbors
			nnIdx,							// nearest neighbors (returned)
			dists,							// distance (returned)
			eps);							// error bound

	// output
	for (int i = 0; i < k; ++i)
	{
		NNPts.push_back(Point_2(dataPts[nnIdx[i]][0], dataPts[nnIdx[i]][1]));
	}
}

#endif // !C_ANN_TREE_WRAPPER_H
