#include "cv_window.h"
#include "../DataColle/voronoi_diagram.h"

int main(int argc, char** argv)
{
	CCVWindow cvWindow(argv[0]);

	// user interface
	while (true)
	{
		// add points
		int x = 0, y = 0;
		int flag;
		// get x, y and flag from the CV window
		flag = cvWindow.getMouseClick(x, y);
		cvWindow.showUpdatedVoronoiDiagram();
	}

	return 0;
}