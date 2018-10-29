#include "mesh_optimizer.h"
#include <CGAL/Polygon_2_algorithms.h> // area


CMeshOptimizer::CMeshOptimizer()
{
}


CMeshOptimizer::~CMeshOptimizer()
{
}

void CMeshOptimizer::SetConstraintPts(std::vector<Point_2> con_pts)
{
	constraint_pts_ = con_pts;
}

void CMeshOptimizer::Update(CVoronoiDiagram *pVD, std::vector<WPoint>& X)
{
	typedef Regular_triangulation::Finite_vertices_iterator RT_Finite_Vertex_Iter;
	typedef Regular_triangulation::Face_handle FaceHandle;
	typedef Regular_triangulation::Vertex_handle VertexHandle;

	Regular_triangulation* rt = pVD->getTriangulation();

	// vertex optimization
	for (RT_Finite_Vertex_Iter viter = rt->finite_vertices_begin(); viter != rt->finite_vertices_end(); ++viter)
	{
		// tranversing the incident faces to the vertex
		VertexHandle vh(viter);
		FaceHandle fh = viter->face();
		FaceHandle fend(fh);
		std::vector<Point_2> face_dual_pos_list;
		std::vector<double> face_area_list;
		do 
		{
			// get index
			int index = 0; // 0, 1, 2
			for (; index < 3; index++)
			{
				if (fh->vertex(index) == vh)
					break;
			}
			// compute
			if (fh->info().size() > 0) 
			{
				face_dual_pos_list.push_back(fh->info().front());
				Point_2 pi = viter->point().point();
				Point_2 pccw = fh->vertex(fh->ccw(index))->point().point();
				Point_2 pcw = fh->vertex(fh->cw(index))->point().point();
				face_area_list.push_back(CGAL::area(pi, pccw, pcw));
			}
			// go to the next face
			FaceHandle fh = fh->neighbor(Regular_triangulation::ccw(index));
		} while (fend != fh);
		// update X with projection constraint
		Point_2 tria_pos_weighted(0, 0);
		double cell_area = 0;
		for (int i = 0; i < face_dual_pos_list.size(); i++)
		{
			tria_pos_weighted += (face_area_list[i]*(face_dual_pos_list[i] - Point_2(0, 0)));
			cell_area += face_area_list[i];
		}
		tria_pos_weighted = Point_2(0, 0) + (tria_pos_weighted - Point_2(0, 0)) / cell_area;
		viter->info() = std::make_pair(tria_pos_weighted, 0.0);
	}
}