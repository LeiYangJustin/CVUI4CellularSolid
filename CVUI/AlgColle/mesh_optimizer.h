#ifndef C_MESH_OPTIMIZER_H
#define C_MESH_OPTIMIZER_H

#include <vector>
#include "../DataColle/voronoi_diagram.h"
#include "algprereq.h"

class ALGCOLLE_CLASS CMeshOptimizer
{
public:
	CMeshOptimizer();
	~CMeshOptimizer();

	void SetConstraintPts(std::vector<Point_2> con_pts);


	// Use the result given by de Goes 2011 HOT to optimize the voronoi decomposition
	//void Update(Regular_triangulation *rt, std::vector<iPoint2> & X);
	void Update(CVoronoiDiagram *vd, std::vector<WPoint> & X);

private:
	std::vector<Point_2> constraint_pts_;
};

#endif // !C_MESH_OPTIMIZER_H



