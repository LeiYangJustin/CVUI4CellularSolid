#ifndef C_RECONSTRUCTOR_H

#include <vector>
#include <set>
#include "datatypedef.h"
#include "cgal_def.h"
#include "voronoi_diagram.h"

class CReconstructor
{
public:
	CReconstructor();
	~CReconstructor();

	void SetReconstructionPts(std::vector<CReconstructionPoint<int>> P_);
	//void Update(Regular_triangulation *rt);
	void Update(CVoronoiDiagram *vd);

private:
	// To-be-reconstructed points
	std::vector< CReconstructionPoint<int> > P_;
	
	void update_edge_by_fitting(std::set<iPoint2>, Edge<int>);
};

#endif // !C_RECONSTRUCTOR_H

