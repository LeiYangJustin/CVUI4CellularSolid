#ifndef C_DATA_TYPE_DEF_H
#define C_DATA_TYPE_DEF_H

#include "prereq.h"

class DATACOLLE_CLASS iPoint2
{
public:
	iPoint2() {};
	iPoint2(int _x, int _y) { x = _x; y = _y; };
	~iPoint2() {};
	int x;
	int y;
};

#endif // !C_DATA_TYPE_DEF_H


