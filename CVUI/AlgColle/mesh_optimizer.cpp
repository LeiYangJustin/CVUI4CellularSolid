#include "mesh_optimizer.h"
#include <CGAL/Polygon_2_algorithms.h> // area
#include <CGAL/centroid.h>
#include "geo_calculator.h"
#include <CGAL/Triangle_2.h>
#include <CGAL/barycenter.h>

CMeshOptimizer::CMeshOptimizer()
{
}


CMeshOptimizer::~CMeshOptimizer()
{
}

void CMeshOptimizer::SetConstraintPts(std::vector<Point_2> con_pts)
{
	for (auto con_pt : con_pts)
	{
		CConstraintPoint cp(con_pt);
		double min_x, min_y, max_x, max_y;
		min_x = con_pt.x() - 20; 
		min_y = con_pt.y() - 20;
		max_x = con_pt.x() + 20; 
		max_y = con_pt.y() + 20;
		std::vector<Point_2> nn_pts;
		for (auto another_con_pt : con_pts)
		{
			if (another_con_pt.x() > min_x && 
				another_con_pt.x() < max_x && 
				another_con_pt.y() > min_y && 
				another_con_pt.y() < max_y)
			{
				nn_pts.push_back(another_con_pt);
			}
		}
		cp.setAllowableDirection(nn_pts);
		constraint_pts_.push_back(cp);

		/*std::cout << con_pt << " has # nn_pts: " << nn_pts.size() 
			<< " and allowable direction: "<< cp.getAllowableDirection() << std::endl;*/
	}
}

void CMeshOptimizer::Update2(CVoronoiDiagram * vd, std::vector<WPoint>& X)
{
	//
	rt_ = vd->getTriangulation();

	//
	update_vertices();

	//
	update_weights();

	//
	std::vector<WPoint> old = X; X.clear();
	for (auto vit = rt_->finite_vertices_begin(); vit != rt_->finite_vertices_end(); ++vit)
	{
		WPoint wp(vit->point().point(), vit->point().weight());
		X.push_back(wp);
	}

	//
	for (int i = 0; i < X.size(); i++)
		std::cout << X[i] << " <----- " << old[i] << std::endl;
}

//void CMeshOptimizer::restore_previous_info()
//{
//	// restore original information at info
//	for (auto vit = rt_->finite_vertices_begin(); vit != rt_->finite_vertices_end(); ++vit)
//		vit->info() = MyWPoint(vit->point().point(), vit->point().weight(), vit->info().is_fixed());
//}

void CMeshOptimizer::update_vertices()
{
	// use c and w to update v
	typedef Regular_triangulation::Finite_vertices_iterator RT_Finite_Vertex_Iter;
	typedef Regular_triangulation::Face_handle FaceHandle;
	typedef Regular_triangulation::Vertex_handle VertexHandle;

	// vertex optimization
	for (RT_Finite_Vertex_Iter viter = rt_->finite_vertices_begin(); viter != rt_->finite_vertices_end(); ++viter)
	{
	//	// restore original information at info
	//	viter->info() = MyWPoint(viter->point().point(), viter->point().weight(), viter->info().is_fixed());

		std::vector<Point_2> face_dual_pos_list;
		std::vector<double> face_area_list;
		Regular_triangulation::Face_circulator fcir = rt_->incident_faces(viter), done(fcir);
		do {
			// compute
			if (!rt_->is_infinite(fcir))
			{
				assert(fcir->info().size() > 0);
				face_dual_pos_list.push_back(fcir->info().front());
				Point_2 p0 = fcir->vertex(0)->point().point();
				Point_2 p1 = fcir->vertex(1)->point().point();
				Point_2 p2 = fcir->vertex(2)->point().point();
				face_area_list.push_back(CGAL::area(p0, p1, p2));
			}
		} while (++fcir != done);

		// compute update
		Vector_2 x_move(0, 0);
		double cell_area = 0;
		for (int i = 0; i < face_dual_pos_list.size(); i++)
		{
			/*
			
			// TAKE WEIGHTS INTO CONSIDERATION

			*/
			x_move += face_area_list[i] * (face_dual_pos_list[i] - Point_2(0, 0));
			cell_area += face_area_list[i];
		}
		x_move = x_move / cell_area;
		x_move -= viter->point().point() - Point_2(0, 0);
		
		// UPDATE X WITH THEIR PROJECTION TO THE SKELETON
		// projecting the motion to allow_dir
		Vector_2 allow_dir = find_allowable_direction(viter->point().point());
		//Vector_2 motion = Point_2(0, 0) + x_new - viter->point().point();
		Vector_2 motion_prj_to_allow_dir = (x_move*allow_dir)*allow_dir;
		Point_2 x_new = find_nearest_constraint_point(viter->point().point() + 0.5*motion_prj_to_allow_dir);

		// handling boundary vertices
		if (!viter->info()) {
			viter->set_point(WPoint(Point_2(x_new.x(), x_new.y()), viter->point().weight()));
		}
	} //end vertex optimization
}

void CMeshOptimizer::update_weights()
{
	// use c and v to update w
	typedef Regular_triangulation::Finite_vertices_iterator RT_Finite_Vertex_Iter;
	typedef Regular_triangulation::Face_handle FaceHandle;
	typedef Regular_triangulation::Vertex_handle VertexHandle;

	// make a vertex indexing
	std::map<int, VertexHandle> map_id_to_vh;
	std::map<VertexHandle, int> map_vh_to_id;
	int num_finite_vertieces = 0;
	for (auto vit = rt_->finite_vertices_begin(); vit != rt_->finite_vertices_end(); ++vit, ++num_finite_vertieces)
	{
		if (map_vh_to_id.find(vit) == map_vh_to_id.end())
		{
			int id = map_vh_to_id.size();
			map_vh_to_id.insert(std::make_pair(vit, id));
			map_id_to_vh.insert(std::make_pair(id, vit));
		}	
	}

	//
	typedef Eigen::Triplet<double> Triplet;
	std::vector<Triplet> triplet_list;
	std::vector<double> dweight_list;
	int num_finite_edges = 0;
	for (auto eit = rt_->finite_edges_begin(); eit != rt_->finite_edges_end(); ++eit, ++num_finite_edges)
	{
		FaceHandle fh_1 = eit->first;
		FaceHandle fh_2 = eit->first->neighbor(eit->second);

		VertexHandle vh_1, vh_2;
		vh_1 = fh_1->vertex(fh_1->cw(eit->second));
		vh_2 = fh_1->vertex(fh_1->ccw(eit->second));

		double dweight = 0.0;
		int num_contributor_to_dweight = 0;
		if (!rt_->is_infinite(fh_1))
		{
			assert(fh_1->info().size() > 0);
			num_contributor_to_dweight++;
			Point_2 pc = fh_1->info().front();

			// dweight = w_2 - w_1
			dweight += (pc - vh_2->point().point()).squared_length() - (pc - vh_1->point().point()).squared_length();
		}
		if (!rt_->is_infinite(fh_2))
		{
			assert(fh_2->info().size() > 0);
			num_contributor_to_dweight++;
			Point_2 pc = fh_2->info().front();
			// dweight = w_2 - w_1
			dweight += (pc - vh_2->point().point()).squared_length() - (pc - vh_1->point().point()).squared_length();
		}
		dweight /= double(num_contributor_to_dweight);

		// add dweight to the vertex indexing
		triplet_list.push_back(Triplet(num_finite_edges, map_vh_to_id[vh_2], 1));
		triplet_list.push_back(Triplet(num_finite_edges, map_vh_to_id[vh_1], -1));
		dweight_list.push_back(dweight);
	}
	// add a constraint // it seems we need to fix a weight in order to solve them all??
	triplet_list.push_back(Triplet(num_finite_edges, 0, 1));
	dweight_list.push_back(0);

	// solve the linear equation constructed from the vertex indexing and the dweight
	// fill A and b
	Eigen::SparseMatrix<double> A(num_finite_edges + 1, num_finite_vertieces);
	A.setFromTriplets(triplet_list.begin(), triplet_list.end());
	Eigen::VectorXd b(num_finite_edges+1), x(num_finite_edges);
	b = Eigen::Map<Eigen::VectorXd>(dweight_list.data(), num_finite_edges + 1, 1);
	// solve Ax = b
	Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> qr_solver;
	qr_solver.compute(A);
	if (qr_solver.info() != Eigen::Success) {
		std::cerr << "decomposition failed" << std::endl;
		return;
	}
	x = qr_solver.solve(b);
	if (qr_solver.info() != Eigen::Success) {
		// solving failed
		std::cerr << "solving failed" << std::endl;
		return;
	}

	// assign
	for (auto vit = rt_->finite_vertices_begin(); vit != rt_->finite_vertices_end(); ++vit)
	{
		int id = map_vh_to_id[vit];
		double w = x(id);
		vit->set_point(WPoint(vit->point().point(), w));
	}
}

Vector_2 CMeshOptimizer::find_allowable_direction(const Point_2 & p)
{
	double min_dist = 10000;
	Vector_2 allow_dir;
	for (auto con_pt : constraint_pts_) {
		double d = (p - con_pt.point_).squared_length();
		if (min_dist > d)
		{
			min_dist = d;
			allow_dir = con_pt.allow_dir_;
		}
	}
	return allow_dir;
}

Point_2 CMeshOptimizer::find_nearest_constraint_point(const Point_2 & p)
{
	double min_dist = 10000;
	Point_2 nearest_con_pt;
	for (auto con_pt : constraint_pts_) {
		double d = (p - con_pt.point_).squared_length();
		if (min_dist > d)
		{
			min_dist = d;
			nearest_con_pt = con_pt.point_;
		}
	}
	return nearest_con_pt;
}
