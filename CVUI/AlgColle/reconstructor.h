#ifndef C_RECONSTRUCTOR_H

#include <vector>
#include <set>
#include "../DataColle/voronoi_diagram.h"
#include "algprereq.h"

class CReconstructionEdge
{
public:
	CReconstructionEdge() {};
	CReconstructionEdge(Point2 p1, Point2 p2) : label_(-1) { p1_ = p1; p2_ = p2; };
	~CReconstructionEdge() {};

	void SetLabel(const int &label) { label_ = label; };
	int GetLabel() const { return label_; };
	bool operator==(const CReconstructionEdge &compare)
	{
		return  compare_points(compare) || swap_compare_points(compare);
	};

	struct edge_less_compare {
		bool operator() (const CReconstructionEdge& lhs, const CReconstructionEdge& rhs) const {
			if (lhs.p1_ != rhs.p1_)
				return lhs.p1_ < rhs.p1_;
			else
				return lhs.p2_ < rhs.p2_;
		}
	};

	Point2 p1_;
	Point2 p2_;

private:
	int label_;

	bool compare_points(const CReconstructionEdge &compare)
	{
		return (p1_ == compare.p1_) && (p2_ == compare.p2_);
	};
	bool swap_compare_points(const CReconstructionEdge &compare)
	{
		return (p1_ == compare.p2_) && (p2_ == compare.p1_);
	};

};


class CReconstructionPoint
{
public:
	CReconstructionPoint() : label_(-1), dist_(0.0), point_(0.0, 0.0) {};
	CReconstructionPoint(const Point2 &_p) : label_(-1), dist_(0.0) { point_ = _p; };
	~CReconstructionPoint() {};
	int GetLabel() { return label_; };
	int GetDist() { return dist_; };
	void SetLabelDist(int label, double dist) { label_ = label; dist_ = dist; };
	double SqrDistToReconEdge(const CReconstructionEdge &e)
	{
		Point2 p1 = e.p1_;
		Point2 p2 = e.p2_;
		Vector2 directVec = e.p2_ - e.p1_;
		Vector2 distVec = (point_ - p1) - ((point_ - p1)*directVec)*directVec / directVec.squared_length();
		return distVec.squared_length();
	}

	struct point_less_compare {
		bool operator() (const CReconstructionPoint& lhs, const CReconstructionPoint& rhs) const {
			return lhs.point_ < rhs.point_;
		}
	};

public:
	Point2 point_;

private:
	int label_;
	double dist_;
};

class ALGCOLLE_CLASS CReconstructor
{
public:
	CReconstructor();
	~CReconstructor();

	void SetReconstructionPts(std::vector<Point2> PointList);
	void Update(CVoronoiDiagram *vd);
	double GetAccuracy();

private:
	double reconstruction_error_;

	// To-be-reconstructed points
	std::vector<CReconstructionPoint> reconPointList_;
	void update_edge_by_fitting(Point2 &fitp1, Point2 &fitp2,
		std::set<CReconstructionPoint, CReconstructionPoint::point_less_compare> rpSet);
};

#endif // !C_RECONSTRUCTOR_H

