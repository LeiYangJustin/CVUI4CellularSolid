#ifndef C_ANN_TREE_WRAPPER_H
#define C_ANN_TREE_WRAPPER_H

#include "algprereq.h"
#include <vector>

#include "../x64_third_party/ANN/ANN/ANN.h"
#include "../DataColle/types.h"

class ANNTreeWrapper
{
public:
	ANNTreeWrapper() : m_eps(0) {};
	~ANNTreeWrapper() { delete m_kdTree; };
	void build_tree(std::vector<Point_2> &dataPts);
	void get_k_nn_from_tree(const Point_2 &in_queryPt, std::vector<Point_2> &NNPts, const int k = 1);
	void get_fr_nn_from_tree(const Point_2 &in_queryPt, std::vector<Point_2> &NNPts, const double sqRadius = 25);

private:
	// Data
	ANNpointArray		m_dataPts;				// data points

	// Query
	double				m_eps;					// error bound

	// Tree
	ANNkd_tree*			m_kdTree;					// search structure
};

#endif // !C_ANN_TREE_WRAPPER_H
