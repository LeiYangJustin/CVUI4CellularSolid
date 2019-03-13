#ifndef _MY_UTILS_
#define _MY_UTILS_

#include <vector>
#include "../LpCVT/LpCVT/common/types.h"

namespace Utils {

	//
	typedef Geex::vec3 vec3;

	//
	static double inner_prdt(vec3 a, vec3 b)
	{
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}

	static double normalize_gradients(std::vector<vec3> &data)
	{
		double norm2 = 0.0;
		for (unsigned i = 0; i < data.size(); ++i)
		{
			//norm2 += data[i] * data[i];
			norm2 += inner_prdt(data[i], data[i]);
		}
		return std::sqrt(norm2);
	};


}

#endif // !_MY_UTILS_
