#ifndef C_VORONOI_DIAGRAM_H
#define C_VORONOI_DIAGRAM_H

#include "prereq.h"
#include "datatypedef.h"
#include <vector>

class DATACOLLE_CLASS CVoronoiDiagram
{
public:
	CVoronoiDiagram();
	~CVoronoiDiagram();

	// add a site with the default weight
	void addSite(WPoint p);
	void getVoronoiSegments(std::vector<std::pair<WPoint, WPoint>> &vseg_list);

private:
	// a VD kernel from CGAL 2D diagram


};

#endif