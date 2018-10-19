#ifndef C_VORONOI_DIAGRAM_H
#define C_VORONOI_DIAGRAM_H

#include "prereq.h"
#include <vector>
//#include "datatypedef.h"
#include "cgal_def.h"

class DATACOLLE_CLASS CVoronoiDiagram
{
public:
	CVoronoiDiagram();
	~CVoronoiDiagram();

	// add a site with the default weight
	void addSite(const BPoint &p);

	// set a list of sites with default weight to triangulation
	void SetSites(const std::vector<BPoint> &pts);

	// pick a site 
	WPoint pickSite(const BPoint &q) const;

	// set the weight to the picked site p
	void setWeightToPickedSite(const BPoint &q, const double &w);

	// get voronoi segments from VD
	void getVoronoiSegments(std::vector<std::pair<WPoint, WPoint>> &vseg_list);

	// return the triangulation
	Regular_triangulation* getTriangulation();

	// update the regular triangulation
	void UpdateTriangulation(const std::vector<WPoint> &wpts);


private:
	// a VD kernel from CGAL 2D diagram
	Regular_triangulation *rt_;

};

#endif