#include "dGBOD12_test.h"
#include "bbnot/scene.h"

bool CTest_dGBOD12::run_test(int nb)
{
	Scene * scene = new Scene;

	// parameters
	int verbose = 0;
	int stepX = 0.0;
	int stepW = 0.0;
	int epsilon = 1.0;
	int frequency = 0;
	int max_iters = 500;

	//
	scene->generate_random_sites(unsigned(nb));

	unsigned iters = scene->optimize_all(stepW, stepX, 500,
		epsilon, max_iters,
		std::cout);

	std::cout << "Finished with iter = " <<  iters << std::endl;

	return true;
}
