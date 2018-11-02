#ifndef C_MESH_OPTIMIZER_H
#define C_MESH_OPTIMIZER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include "algprereq.h"

#include "../DataColle/voronoi_diagram.h"

// ???
#include "../DataColle/reconstruction_data.h"

class ALGCOLLE_CLASS CMeshOptimizer
{
public:
	CMeshOptimizer();
	~CMeshOptimizer();

	void SetConstraintPts(std::vector<Point_2> con_pts);


	// Use the result given by de Goes 2011 HOT to optimize the voronoi decomposition
	//void Update(Regular_triangulation *rt, std::vector<iPoint2> & X);
	//void Update(CVoronoiDiagram *vd, std::vector<WPoint> & X);

	void Update2(CVoronoiDiagram *vd, std::vector<WPoint> & X);

private:
	std::vector<CConstraintPoint> constraint_pts_;

	Regular_triangulation * rt_;

	//void restore_previous_info();
	void update_vertices();
	void update_weights();

	Vector_2 find_allowable_direction(const Point_2 &p);
	Point_2 find_nearest_constraint_point(const Point_2 &p);
};

#endif // !C_MESH_OPTIMIZER_H



