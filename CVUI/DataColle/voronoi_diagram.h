#ifndef C_VORONOI_DIAGRAM_H
#define C_VORONOI_DIAGRAM_H

#include "prereq.h"
#include <vector>
#include "datatypedef.h"
#include "cgal_def.h"

class DATACOLLE_CLASS CVoronoiDiagram
{
public:
	CVoronoiDiagram();
	~CVoronoiDiagram();

	// add a site with the default weight
	void addSite(const iPoint2 &p);

	// pick a site 
	WPoint pickSite(const iPoint2 &q);

	// set the weight to the picked site p
	void setWeightToPickedSite(const BPoint &q, const double &w);

	// get voronoi segments from VD
	void getVoronoiSegments(std::vector<std::pair<iPoint2, iPoint2>> &vseg_list);

private:
	// a VD kernel from CGAL 2D diagram
	Regular_triangulation *rt_;

};

#endif