// this is adapted from src code of de Goes et al. (2012)
// Blue Noise through Optimal Transport
// from http://fernandodegoes.org/

#ifndef _RANDOM_
#define _RANDOM_ 1

inline
double random_double(const double min, const double max)
{
	double range = max - min;
	return min + (double(rand()) / double(RAND_MAX)) * range;
}

inline
int random_int(const int min, const int max)
{
    int range = max - min;
    return min + int((double(rand())/double(RAND_MAX)) * range);
}

template <class Vector>
Vector random_vec(const double scale)
{
    double dx = random_double(-scale, scale);
    double dy = random_double(-scale, scale);
    return Vector(dx, dy);
}

#endif
