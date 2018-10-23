#ifndef C_RECONSTRUCTOR_H
#define C_RECONSTRUCTOR_H

#include <vector>
#include <set>
#include "../DataColle/voronoi_diagram.h"
#include "algprereq.h"

//class CReconstructionEdge
//{
//public:
//	CReconstructionEdge():label_(-1) {};
//	CReconstructionEdge(Point_2 p1, Point_2 p2) : label_(-1) { p1_ = p1; p2_ = p2; };
//	~CReconstructionEdge() {};
//
//	void SetLabel(const int &label) { label_ = label; };
//	int GetLabel() const { return label_; };
//	bool operator==(const CReconstructionEdge &compare)
//	{
//		return  compare_points(compare) || swap_compare_points(compare);
//	};
//
//	struct edge_less_compare {
//		bool operator() (const CReconstructionEdge& lhs, const CReconstructionEdge& rhs) const {
//			if (lhs.p1_ != rhs.p1_)
//				return lhs.p1_ < rhs.p1_;
//			else
//				return lhs.p2_ < rhs.p2_;
//		}
//	};
//
//	Point_2 p1_;
//	Point_2 p2_;
//
//private:
//	int label_;
//
//	bool compare_points(const CReconstructionEdge &compare)
//	{
//		return (p1_ == compare.p1_) && (p2_ == compare.p2_);
//	};
//	bool swap_compare_points(const CReconstructionEdge &compare)
//	{
//		return (p1_ == compare.p2_) && (p2_ == compare.p1_);
//	};
//
//};


class CReconstructionPoint
{
public:
	CReconstructionPoint() : label_(-1), dist_(100.0), point_(0.0, 0.0) {};
	CReconstructionPoint(const Point_2 &_p) : label_(-1), dist_(100.0) { point_ = _p; };
	~CReconstructionPoint() {};
	int GetLabel() { return label_; };
	int GetDist() { return dist_; };
	void SetLabelDist(int label, double dist) { label_ = label; dist_ = dist; };
	//double SqrDistToReconEdge(const CReconstructionEdge &e)
	//{
	//	Point_2 p1 = e.p1_;
	//	Point_2 p2 = e.p2_;
	//	Vector_2 directVec = e.p2_ - e.p1_;
	//	Vector_2 distVec = (point_ - p1) - ((point_ - p1)*directVec)*directVec / directVec.squared_length();
	//	return distVec.squared_length();
	//}
	double SqrDistToSegment(const Segment_2 &e)
	{
		Vector_2 a = point_ - e.source();
		Vector_2 v = e.target() - e.source();
		Vector_2 proj_a_v = v*(a*v / (v*v));
		double t = sqrt(proj_a_v.squared_length()) / sqrt(v.squared_length());
		double s = a*v;
		if (s < 0)
			t = -t;
		if (t >= 0 && t <= 1)
		{
			Vector_2 d = a - proj_a_v;
			return sqrt(d.squared_length());
		}
		else if (t < 0) {
			return sqrt(a.squared_length());
		}
		else {
			return sqrt((point_ - e.target()).squared_length());
		}
	}

	struct point_less_compare {
		bool operator() (const CReconstructionPoint& lhs, const CReconstructionPoint& rhs) const {
			return lhs.point_ < rhs.point_;
		}
	};

public:
	Point_2 point_;

private:
	int label_;
	double dist_;
};

class ALGCOLLE_CLASS CReconstructor
{
public:
	CReconstructor();
	~CReconstructor();

	void SetBBox(int w, int h);
	void SetReconstructionPts(std::vector<Point_2> PointList);
	void GetReconstructionPts(std::vector<Point_2> &RPts, std::vector<int> &RPt_labels);
	void Update(CVoronoiDiagram *vd);
	double GetAccuracy();

private:
	int width_;
	int height_;
	double reconstruction_error_;

	// To-be-reconstructed points
	std::vector<CReconstructionPoint> reconPointList_;
	void update_edge_by_fitting(Point_2 &fitp1, Point_2 &fitp2,
		std::list<CReconstructionPoint> rp_list);

	void update_edge_by_fitting_new(Point_2 &fitp1, Point_2 &fitp2,
		std::list<CReconstructionPoint> rp_list);
};

#endif // !C_RECONSTRUCTOR_H

