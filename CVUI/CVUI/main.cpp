#include "cv_window.h"
#include "../DataColle/datatypedef.h"
#include "../DataColle/voronoi_diagram.h"

int main(int argc, char** argv)
{
	CCVWindow cvWindow(argv[0]);

	// user interface
	while (true)
	{
		// add points
		int x =0 , y = 0;
		// get x and y from the CV window
		cvWindow.getMouseClickForAddingVSite(x, y);

		// update the new VD in the window
		cvWindow.showVoronoiDiagram();
	}

	return 0;
}