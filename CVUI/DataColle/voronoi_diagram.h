#ifndef C_VORONOI_DIAGRAM_H
#define C_VORONOI_DIAGRAM_H

#include "dataprereq.h"
#include <vector>
#include <iostream>
#include "cgal_def.h"



struct Cropped_voronoi_from_delaunay {
	std::list<Segment_2> m_cropped_vd;
	Iso_rectangle_2 m_bbox;

	Cropped_voronoi_from_delaunay(const Iso_rectangle_2& bbox) :m_bbox(bbox) {}

	template <class RSL>
	void crop_and_extract_segment(const RSL& rsl) {
		CGAL::Object obj = CGAL::intersection(rsl, m_bbox);
		const Segment_2* s = CGAL::object_cast<Segment_2>(&obj);
		if (s) m_cropped_vd.push_back(*s);
	}

	void operator<<(const Ray_2& ray) { crop_and_extract_segment(ray); }
	void operator<<(const Line_2& line) { crop_and_extract_segment(line); }
	void operator<<(const Segment_2& seg) { crop_and_extract_segment(seg); }
};

class DATACOLLE_CLASS CVoronoiDiagram
{
public:
	CVoronoiDiagram();
	~CVoronoiDiagram();

	// return the triangulation
	Regular_triangulation* getTriangulation();

	// update the regular triangulation
	void UpdateTriangulation(const std::vector<WPoint> &wpts);

	void GetCroppedVoronoiSegments(Iso_rectangle_2 bbox, std::vector<Segment_2> &voronoi_segments);
	void GetFittedSegments(std::vector<Segment_2> &fitting_segments);
	void GetDualEdges(std::vector<Segment_2> &dual_edges);
	void GetFittingBasePtsMap(std::map<Point_2, std::vector<Point_2>> &fitting_base_pts_map);
private:
	// a VD kernel from CGAL 2D diagram
	Regular_triangulation* rt_;

	// storage cropped segments
	std::vector<Segment_2> voronoi_segments_;

	void clearMembers() {
		rt_->clear();
		voronoi_segments_.clear();
	};
};

#endif