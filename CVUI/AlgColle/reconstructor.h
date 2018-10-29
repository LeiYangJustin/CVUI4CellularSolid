#ifndef C_RECONSTRUCTOR_H
#define C_RECONSTRUCTOR_H

#include <vector>
#include <set>
#include "../DataColle/voronoi_diagram.h"
#include "algprereq.h"
#include "../DataColle/reconstruction_data.h"

class ALGCOLLE_CLASS CReconstructor
{
public:
	CReconstructor();
	~CReconstructor();

	void SetBBox(int w, int h);
	void SetReconstructionPts(std::vector<Point_2> PointList);
	void GetReconstructionPts(std::vector<Point_2> &RPts, std::vector<int> &RPt_labels);
	void Update(CVoronoiDiagram *vd);
	void Update2(CVoronoiDiagram *vd);
	double GetReconstructionError();

private:
	int width_;
	int height_;
	double reconstruction_error_;

	// To-be-reconstructed points
	std::vector<CReconstructionPoint> reconPointList_;

	void update_edge_by_fitting(Point_2 &fitp1, Point_2 &fitp2,
		std::list<CReconstructionPoint> rp_list);

	double update_edge_by_fitting_new(Point_2 &fitp1, Point_2 &fitp2,
		std::list<CReconstructionPoint> rp_list);

	//void assign_label_to_reconstruction_pts(CVoronoiDiagram *vd,
	//	std::map<int, std::list<CReconstructionPoint>> &map_lid_to_rplist);

	void fit_line_segment_to_labeled_pts(CVoronoiDiagram *pVD,
		std::vector<Segment_2> voronoi_segments,
		std::map<int, std::list<CReconstructionPoint>> map_lid_to_rplist);

};

#endif // !C_RECONSTRUCTOR_H

