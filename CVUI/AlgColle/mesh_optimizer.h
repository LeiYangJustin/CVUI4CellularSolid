#ifndef C_MESH_OPTIMIZER_H

#include <vector>
#include "../DataColle/voronoi_diagram.h"
#include "algprereq.h"

class ALGCOLLE_CLASS CMeshOptimizer
{
public:
	CMeshOptimizer();
	~CMeshOptimizer();

	// Use the result given by de Goes 2011 HOT to optimize the voronoi decomposition
	//void Update(Regular_triangulation *rt, std::vector<iPoint2> & X);
	void Update(CVoronoiDiagram *vd, std::vector<WPoint> & X);
};

#endif // !C_MESH_OPTIMIZER_H



