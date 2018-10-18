#ifndef C_DATA_TYPE_DEF_H
#define C_DATA_TYPE_DEF_H

#include "prereq.h"

// base
template <class DataType> class Point
{
public:
	DataType x;
	DataType y;

public:
	Point() {};
	virtual ~Point() {};

	/// Overloading operators
	// overload operator=
	void operator=(Point &p)
	{
		x = p.x;
		y = p.y;
	};
	// overload operator==
	bool operator==(Point &p)
	{
		if (abs(x - p.x)<0.000000001 && abs(y - p.y)< 0.000000001)
			return true;
		else
			return false;
	}
	// overload operator+ for computing the sum
	Point operator+(Point &add)
	{
		Point p;
		p.x = x + add.x;
		p.y = y + add.y;
		return p;
	};
	// overload operator-
	Point operator-(Point &minus)
	{
		Point p;
		p.x = x - minus.x;
		p.y = y - minus.y;
		return p;
	};
	// overload operator/ for computing the average
	Point operator/(const float divider)
	{
		Point p;
		p.x = x/divider;
		p.y = y/divider;
		return p;
	};
	double norm()
	{
		return sqrt(x*x + y*y);
	}
};

// define Edge as a two-point object
template <class DataType> class Edge
{
public:
	Edge() {};
	Edge(Point<DataType> p1, Point<DataType> p2)
	{
		first = p1;
		second = p2;
	}
	Edge(Edge &e)
	{
		first = e.first;
		second = e.second;
	}
	~Edge() {};

public:
	Point<DataType> first;
	Point<DataType> second;
};

// this is for OpenCV visualization
class DATACOLLE_CLASS iPoint2: public Point<int>
{
public:
	iPoint2() {};
	iPoint2(int _x, int _y)
	{
		x = _x; y = _y; 
	};
	~iPoint2() {};
	int x;
	int y;
};

// this is a templated class to store the to-be-reconstructed point
template <class DataType> class  DATACOLLE_CLASS CReconstructionPoint : public Point<DataType>
{
private:
	std::pair<int, double>  label_dist_;

public:
	CReconstructionPoint()
	{
		x = 0;
		y = 0;
		label_dist_ = std::pair<int, double>(-1, 0);
	};

	CReconstructionPoint(int _x, int _y)
	{
		x = _x;
		y = _y;
		label_dist_ = std::pair<int, double>(-1, 0);
	};

	CReconstructionPoint(CReconstructionPoint &p)
	{
		x = p.x;
		y = p.y;
		if (p.GetLabelAndMinDist().first != -1)
			label_dist_ = p.GetLabelAndMinDist();
	};

	~CReconstructionPoint() {};

	/// Label and min distance
	// set label and min_dist
	void SetLabelAndMinDist(const int &l, const double &d)
	{
		if (label_dist_.first == -1) {
			label_dist_ = std::pair<int, double>(l, d);
		}
		else if (label_dist_.second > d) {
			label_dist_ = std::pair<int, double>(l, d);
		}
	};
	// get label and min dist
	std::pair<int, double> GetLabelAndMinDist() const
	{
		return label_dist_;
	};
	// get label
	int GetLabel() const
	{
		return label_dist_.first;
	};
	// get minimal distance
	double GetMinDist() const
	{
		if (label_dist_.first != -1)
			return label_dist_.second;
		else
			return 0.0;
	};

	/// Point to X distance
	// point to edge distance
	double DistToEdge(Edge<int> e)
	{
		double d = 0.0;
		/* compute distance from this to Edge e */
		return d;
	};
	double DistToPoint(CReconstructionPoint e)
	{
		double d = 0.0;
		/* compute distance from this to Edge e */
		return d;
	};
};

//// this is a templated class to describe the vertices and sites of the voronoi diagram
//template <class DataType> class  DATACOLLE_CLASS CVDPoint : public Point<DataType>
//{
//private:
//	KPoint kp;
//
//public:
//	CVDPoint() {};
//	CVDPoint(KPoint &_kp)
//	{
//		kp = _kp;
//		x = kp.x;
//		y = kp.y;
//	};
//	~CVDPoint() {};
//
//	/* see what we need */
//};

#endif // !C_DATA_TYPE_DEF_H


