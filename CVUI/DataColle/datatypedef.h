#ifndef C_DATA_TYPE_DEF_H
#define C_DATA_TYPE_DEF_H
#include "prereq.h"

class DATACOLLE_CLASS WPoint
{
public:
	WPoint(int _x, int _y):w_(0.0) 
	{
		x_ = _x; 
		y_ = _y; 
	};
	WPoint(int _x, int _y, double _w) 
	{ 
		x_ = _x; 
		y_ = _y; 
		w_ = _w; 
	};
	~WPoint() {};

	// get data
	int x() const { return x_; };
	int y() const { return y_; };
	double weight() const { return w_; };

	// set data
	void setX(const int _x) { x_ = _x; };
	void setY(const int _y) { y_ = _y; };
	void setWeight(const double _w) { w_ = _w; };

private:
	int x_;
	int y_;
	double w_;
};

#endif // !C_DATA_TYPE_DEF_H


