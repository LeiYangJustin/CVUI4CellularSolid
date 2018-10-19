#include "reconstructor.h"
#include <map>


CReconstructor::CReconstructor() :
	reconstruction_error_(0.0)
{
}


CReconstructor::~CReconstructor()
{
}

void CReconstructor::SetReconstructionPts(std::vector<Point2> P)
{
	for (int i = 0; i < P.size(); i++)
	{
		CReconstructionPoint rp(P[i]);
		reconPointList_.push_back(rp);
	}
}

void CReconstructor::Update(CVoronoiDiagram *vd)
{
	typedef Regular_triangulation::Finite_edges_iterator RT_Finite_Edge_Iter;
	typedef Regular_triangulation::Face_handle FaceHandle;
	typedef Regular_triangulation::Vertex_handle VertexHandle;
	
	// update rt in vd
	Regular_triangulation *rt = vd->getTriangulation();

	// extract edges by filtering mirrored halfedges
	std::map<CReconstructionEdge, std::set<CReconstructionPoint>> map_edges_to_rpset;
	for (RT_Finite_Edge_Iter eit = rt->finite_edges_begin(); eit != rt->finite_edges_end(); ++eit)
	{
		// point to edge computation that returns a distance
		FaceHandle fh = eit->first;
		int vind = eit->second;
		Point2 pstart = fh->vertex(vind + 1)->point().point();
		Point2 pend = fh->vertex(vind + 2)->point().point();
		CReconstructionEdge tmpEdge(pstart, pend);

		// check availabitility
		bool b_availability = true;
		for (auto miter = map_edges_to_rpset.begin(); miter != map_edges_to_rpset.end(); ++miter)
		{
			if (tmpEdge == miter->first)
			{
				b_availability = false;
				break;
			}
		}
		if (b_availability) {
			tmpEdge.SetLabel(map_edges_to_rpset.size());
			std::set<CReconstructionPoint> emptySet;
			map_edges_to_rpset.insert(std::make_pair(tmpEdge, emptySet));
		}
	}

	// assign label and distance to each reconstruction point
	// make point sets according to the labels
	// compute reconstruction_error_
	reconstruction_error_ = 0.0;
	for (auto rp_iter = reconPointList_.begin(); rp_iter != reconPointList_.end(); ++rp_iter)
	{
		CReconstructionEdge tmpEdge;
		// assign label and distance to each reconstruction point
		for (auto miter = map_edges_to_rpset.begin(); miter != map_edges_to_rpset.end(); ++miter)
		{
			double dist = rp_iter->SqrDistToReconEdge(miter->first);
			if (rp_iter->GetDist() > dist)
			{
				rp_iter->SetLabelDist(miter->first.GetLabel(), dist);
				tmpEdge = miter->first;
			}
		}

		// add the point to the point set
		if (map_edges_to_rpset.find(tmpEdge) != map_edges_to_rpset.end())
		{
			map_edges_to_rpset[tmpEdge].insert(*rp_iter);
		}
		else
		{
			// Oops! something wrong
			std::cerr << "Oops! something wrong" << std::endl;
		}

		// compute reconstruction_error_
		reconstruction_error_ += rp_iter->GetDist();
	}

	// fit a line segment to each of the point sets
	for (auto miter = map_edges_to_rpset.begin(); miter != map_edges_to_rpset.end(); ++miter)
	{
		Point2 fitp1, fitp2;
		update_edge_by_fitting(miter->second, fitp1, fitp2);
		VertexHandle vh1 = rt->nearest_power_vertex(miter->first.p1_);
		vh1->info().push_back(fitp1);
		VertexHandle vh2 = rt->nearest_power_vertex(miter->first.p2_);
		vh1->info().push_back(fitp2);
	}
}

double CReconstructor::GetAccuracy()
{
	return reconstruction_error_;
}

void CReconstructor::update_edge_by_fitting(std::set<CReconstructionPoint> rpSet, Point2 fitp1, Point2 fitp2)
{
	// fitting a line in the first place;
	// PCA

	// finding the two ends
	// Projection

}
