#include "ann_tree_wrapper.h"

void ANNTreeWrapper::build_tree(std::vector<Point_2> &in_dataPts)
{
	int nPts = in_dataPts.size();
	int dim = 2;
	m_dataPts = annAllocPts(nPts, dim);			// allocate data points

	for (int i = 0; i < in_dataPts.size(); ++i)
	{
		m_dataPts[i][0] = in_dataPts[i].x();
		m_dataPts[i][1] = in_dataPts[i].y();
	}

	m_kdTree = new ANNkd_tree(		// build search structure
		m_dataPts,					// the data points
		nPts,						// number of points
		dim);						// dimension of space
}

void ANNTreeWrapper::get_k_nn_from_tree(const Point_2 &in_queryPt, std::vector<Point_2> &NNPts, int k)
{
	auto queryPt = annAllocPt(2);					// allocate query point
	auto nnIdx = new ANNidx[k];						// allocate near neigh indices
	auto dists = new ANNdist[k];					// allocate near neighbor dists

	queryPt[0] = in_queryPt.x();
	queryPt[1] = in_queryPt.y();

	m_kdTree->annkSearch(						// search
		queryPt,						// query point
		k,								// number of near neighbors
		nnIdx,							// nearest neighbors (returned)
		dists,							// distance (returned)
		m_eps);							// error bound
	
	for (int i = 0; i < k; ++i)			// output
	{
		NNPts.push_back(Point_2(m_dataPts[nnIdx[i]][0], m_dataPts[nnIdx[i]][1]));
	}

	delete nnIdx;
	delete dists;
}

void ANNTreeWrapper::get_fr_nn_from_tree(const Point_2 & in_queryPt, std::vector<Point_2>& NNPts, const double sqRadius)
{
	auto queryPt = annAllocPt(2);					// allocate query point
	int k = 100;
	auto nnIdx = new ANNidx[k];						// allocate near neigh indices
	auto dists = new ANNdist[k];					// allocate near neighbor dists

	queryPt[0] = in_queryPt.x();
	queryPt[1] = in_queryPt.y();

	m_kdTree->annkFRSearch(			// approx fixed-radius kNN search
		queryPt,				// query point
		sqRadius,				// squared radius
		k,						// number of near neighbors to return
		nnIdx,					// nearest neighbor array (modified)
		dists,					// dist to near neighbors (modified)
		m_eps);					// error bound

	//size_t sz = nnIdx

	for (int i = 0; i < k; ++i)			// output
	{
		if (nnIdx[i] == ANN_NULL_IDX)
			break;
		//std::cout << m_dataPts[nnIdx[i]][0] << ", " << m_dataPts[nnIdx[i]][1] << std::endl;
		NNPts.push_back(Point_2(m_dataPts[nnIdx[i]][0], m_dataPts[nnIdx[i]][1]));
	}

	delete nnIdx;
	delete dists;
}
