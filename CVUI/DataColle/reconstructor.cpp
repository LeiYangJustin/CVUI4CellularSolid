#include "reconstructor.h"



CReconstructor::CReconstructor()
{
}


CReconstructor::~CReconstructor()
{
}

void CReconstructor::SetReconstructionPts(std::vector<CReconstructionPoint<int>> P)
{
	P_ = P;
}

void CReconstructor::Update(CVoronoiDiagram *vd)
{
	/*
	this seems not correct; check later
	*/
	Regular_triangulation *rt = new Regular_triangulation;
	vd->getTriangulation(rt);

	typedef Regular_triangulation::Finite_edges_iterator RT_Finite_Edge_Iter;
	for (RT_Finite_Edge_Iter eit = rt->finite_edges_begin(); eit != rt->finite_edges_end(); ++eit)
	{
		std::set<iPoint2> pset;
		for (int i = 0; i < P_.size(); i++)
		{
			// point to edge computation that returns a distance

			// if condition is met
			// assign label-mindist to the pt and add to pset 

		}

		// fit a line segment to each of the point sets
		Edge<int> edge;
		update_edge_by_fitting(pset, edge);

		// update the eit data with edge

	}

	// finally we compute the averaged position of the updated vertices by averaging all edge ends related to the vertices
	// update rt
	// check vd

}

void CReconstructor::update_edge_by_fitting(std::set<iPoint2>, Edge<int>)
{
	// do something
}
