#ifndef C_VORONOI_DIAGRAM_H
#define C_VORONOI_DIAGRAM_H

#include "dataprereq.h"
#include <vector>
#include "cgal_def.h"

class DATACOLLE_CLASS CVoronoiDiagram
{
public:
	CVoronoiDiagram();
	~CVoronoiDiagram();

	// return the triangulation
	Regular_triangulation* getTriangulation();

	// update the regular triangulation
	void UpdateTriangulation(const std::vector<WPoint> &wpts);

	void GetCroppedVoronoiSegments(std::vector<std::pair<Point2, Point2>> &voronoi_segments);

private:
	// a VD kernel from CGAL 2D diagram
	Regular_triangulation *rt_;





	//// add a site with the default weight
	//void AddSite(const BPoint &p);

	//// set a list of sites with default weight to triangulation
	//void SetSites(const std::vector<BPoint> &pts);

	//// pick a site 
	//WPoint PickSite(const BPoint &q) const;

	//// set the weight to the picked site p
	//void SetWeightToPickedSite(const BPoint &q, const double &w);

	//// get voronoi segments from VD
	//void GetVoronoiSegments(std::vector<std::pair<WPoint, WPoint>> &vseg_list);

};

#endif