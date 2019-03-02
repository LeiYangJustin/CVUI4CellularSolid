// this is adapted from src code of de Goes et al. (2012)
// Blue Noise through Optimal Transport
// from http://fernandodegoes.org/

#include <map>
#include "scene.h"
//#include "timer.h"
#include "pw_line_search.h"

typedef CLSWeights<Scene, FT> LSWeights;
typedef CLSPositions<Scene, Point_2, Vector_2> LSPositions;

FT Scene::optimize_positions_via_lloyd(bool update)
{
    std::vector<Point_2> points;
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;        
		Point_2 ci = m_rt.compute_centroid(vi);
		points.push_back(ci);

		//if (m_domain.is_outside(vi->get_position()))
		//	points.push_back(vi->get_position());
		//else{
		//	Point_2 ci = m_rt.compute_centroid(vi);
		//	points.push_back(ci);
		//}
    }
    update_positions(points);
    if (update) update_triangulation();

	//update_positions(points, false, true);
	//if (update) update_triangulation();

    std::vector<Vector_2> gradient;
    compute_position_gradient(gradient);
    return compute_norm(gradient);
}

FT Scene::optimize_positions_via_gradient_ascent(FT timestep, bool update)
{
    std::vector<Point_2> points;
    collect_visible_points(points);
    
    std::vector<Vector_2> gradient;
    compute_position_gradient(gradient);
    
    if (timestep <= 0.0)
    {
        double mean_capacity = compute_mean(m_capacities);
        double max_alpha = 1.0 / mean_capacity;        
        LSPositions line_search(this, 100, max_alpha);
        
        //if (m_timer_on) { Timer::start_timer(m_timer, COLOR_RED, "XLineSearch"); std::cout << std::endl; }
        FT step = line_search.run_bt(points, gradient);
        //if (m_timer_on) Timer::stop_timer(m_timer, COLOR_RED);
        
        gradient.clear();
        compute_position_gradient(gradient);
        return compute_norm(gradient);
    }

    for (unsigned i = 0; i < points.size(); ++i)
    {
        Point_2  pi = points[i];
        Vector_2 gi = gradient[i];
        points[i] = pi + timestep*gi;
    }

    update_positions(points);
    if (update) update_triangulation();
    
    gradient.clear();
    compute_position_gradient(gradient);
    return compute_norm(gradient);
}

FT Scene::optimize_weights_via_gradient_descent(FT timestep, bool update)
{
    std::vector<FT> gradient;
    compute_weight_gradient(gradient, -1.0);
    
    std::vector<FT> weights;
    collect_visible_weights(weights);
    
    if (timestep <= 0.0)
    {
        LSWeights line_search(this, 20, 2.0);
        FT step = line_search.run_bt(weights, gradient);
    } else {
        for (unsigned i = 0; i < weights.size(); ++i)
        {
            FT wi = weights[i];
            FT gi = gradient[i];
            weights[i] = wi + timestep*gi;
        }
        update_weights(weights);
        if (update) update_triangulation();
    }

    compute_weight_gradient(gradient);
    return compute_norm(gradient);
}

FT Scene::optimize_weights_via_newton(FT timestep, bool update)
{
    std::vector<FT> gradient;
    compute_weight_gradient(gradient, -1.0);
    
    std::vector<FT> direction;
    bool ok = solve_newton_step(gradient, direction);
    if (!ok) return 0.0;
    
    std::vector<FT> weights;
    collect_visible_weights(weights);
    
    if (timestep <= 0.0)
    {
        LSWeights line_search(this, 20, 2.0);
        FT step = line_search.run_bt(weights, direction);
    } else {
        for (unsigned i = 0; i < weights.size(); ++i)
        {
            FT wi = weights[i];
            FT gi = direction[i];
            weights[i] = wi + timestep*gi;
        }    
        update_weights(weights);
        if (update) update_triangulation();
    }
    
    compute_weight_gradient(gradient);
    return compute_norm(gradient);
}

bool Scene::solve_newton_step(const std::vector<FT>& b, std::vector<FT>& x)
{
    //if (m_timer_on) Timer::start_timer(m_timer, COLOR_BLUE, "LinearSolver");

    unsigned nb = 0;
    std::map<unsigned, unsigned> indices;
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;        
        indices[vi->get_index()] = nb++;
    }
    
    Eigen::SparseMatrix<double> L(nb, nb);
    build_laplacian(0.5, indices, L);
    
    bool ok = solve_linear_system(L, x, b);
    if (!ok) 
    {
        //std::cout << red << "linear solver failed" << white << std::endl;
		std::cout << "linear solver failed" << std::endl;
        return false;
    }

    //if (m_timer_on) Timer::stop_timer(m_timer, COLOR_BLUE);
    return true;
}

void Scene::build_laplacian(const FT scale,
                            const std::map<unsigned, unsigned>& indices,
							Eigen::SparseMatrix<double>& A) const
{
    //unsigned nb = A.rows();

	std::vector<Eigen::Triplet<double>> triplet_list;

    for (unsigned k = 0; k < m_vertices.size(); ++k)
    {
        Vertex_handle vi = m_vertices[k];
        if (vi->is_hidden()) continue;
        unsigned i = indices.find(vi->get_index())->second;
        
        double diagi = 0.0;
        //SparseArray rowi(nb);        
        Edge_circulator ecirc = m_rt.incident_edges(vi);
        Edge_circulator eend  = ecirc;
        CGAL_For_all(ecirc, eend)
        {
            Edge edge = *ecirc;
            if (!m_rt.is_inside(edge)) continue;
            
            Vertex_handle vj = m_rt.get_source(edge);
            if (vj == vi) vj = m_rt.get_target(edge);
            
            unsigned j = vj->get_index();
            j = indices.find(j)->second;
            
            double coef = scale * m_rt.get_ratio(edge);
            if (std::abs(coef) < MY_EPS) continue;
            
			triplet_list.push_back(Eigen::Triplet<double>(i, j, -coef));
            //rowi.setValue(j, -coef);
            diagi += coef;
        }
		triplet_list.push_back(Eigen::Triplet<double>(i, i, diagi));

        //rowi.setValue(i, diagi);
        //A.setRow(i, rowi);
    }
	A.setFromTriplets(triplet_list.begin(), triplet_list.end());
}

bool Scene::solve_linear_system(const Eigen::SparseMatrix<double>& A,
                                std::vector<double>& x,
								const std::vector<double>& b) const
{
	//// original code
    //SuiteSparseQRFactorizer solver;
    //bool ok = solver.factorize(A);
    //if (!ok) return false;
    //
    //ok = solver.solve(b, x);
    //return ok;

	auto bvector = b;
	Eigen::VectorXd bb = Eigen::Map<Eigen::VectorXd>(bvector.data(), bvector.size());

	// fill b
	// solve Ax = b
	Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> QR_solver;
	QR_solver.compute(A);
	if (QR_solver.info() != Eigen::Success) {
		// decomposition failed
		std::cout << "decomposition failed" << std::endl;
		return false;
	}
	Eigen::VectorXd xx = QR_solver.solve(bb);
	if (QR_solver.info() != Eigen::Success) {
		// solving failed
		std::cout << "solving failed" << std::endl;
		return false;
	}

	// output x
	x = std::vector<double>(xx.data(), xx.data() + xx.rows() * xx.cols());
	return true;
}

unsigned Scene::optimize_weights_via_gradient_descent_until_converge(FT timestep, 
                                                                     FT threshold,
                                                                     unsigned update,
                                                                     unsigned max_iters)
{
    for (unsigned i = 0; i < max_iters; ++i)
    {
        bool flag = (update == 0 || (i+1) % update == 0);
        FT norm = optimize_weights_via_gradient_descent(timestep, flag);
        if (norm < threshold) return i;
    }
    return max_iters;
}

unsigned Scene::optimize_weights_via_newton_until_converge(FT timestep, 
                                                           FT threshold,
                                                           unsigned update,
                                                           unsigned max_iters)
{
    for (unsigned i = 0; i < max_iters; ++i)
    {
        bool flag = (update == 0 || (i+1) % update == 0);
        FT norm = optimize_weights_via_newton(timestep, flag);
        if (norm < threshold) return i;
    }
    return max_iters;
}

unsigned Scene::optimize_all(FT wstep, FT xstep, 
                             unsigned max_newton_iters,
                             FT epsilon, 
                             unsigned max_iters,
                             std::ostream& out)
{
    bool global_connectivity = m_fixed_connectivity;
    FT norm = optimize_positions_via_lloyd(true);    
    unsigned nb0 = count_visible_sites();
    
    FT xthreshold = compute_position_threshold(epsilon);
    FT wthreshold = compute_weight_threshold(epsilon);

    out << "NbSites: " << nb0 << std::endl;
    out << "Threshold: " << xthreshold << " ; " << wthreshold << std::endl;

	std::string fname = "synthss_";
	std::string fname_increment = fname + "0";
	std::string affix = ".mesh";
	write_updated_triangulation(fname_increment +affix, m_rt);

    m_fixed_connectivity = false;
    FT coarse_xthreshold = 2.0*xthreshold;
    FT coarse_wthreshold = 2.0*wthreshold;

    unsigned iters = 0;
    while (iters < max_iters)
    {
        iters++;
        reset_weights();        
        optimize_weights_via_newton_until_converge(wstep, coarse_wthreshold, 0, max_newton_iters);        
        norm = optimize_positions_via_lloyd(true);

		//
		fname_increment = fname + std::to_string(iters);
		write_updated_triangulation(fname_increment + affix, m_rt);
		//

        std::cout << "Norm: " << norm << std::endl;
        if (norm <= coarse_xthreshold) break;
    }
    
    std::cout << "Partial: " << iters << " iters" << std::endl;
    m_fixed_connectivity = global_connectivity;
    if (iters == max_iters) return iters;
    
    m_fixed_connectivity = false;
    FT fine_xthreshold = xthreshold;
    FT fine_wthreshold = wthreshold;
    
    while (iters < max_iters)
    {
        iters++;        
        unsigned nb1 = count_visible_sites();
        if (nb1 != nb0) reset_weights();
        optimize_weights_via_newton_until_converge(wstep, fine_wthreshold, 0, max_newton_iters);
        norm = optimize_positions_via_gradient_ascent(xstep, true);

		//
		fname_increment = fname + std::to_string(iters);
		write_updated_triangulation(fname_increment + affix, m_rt);
		//

        std::cout << "Norm: " << norm << std::endl;
        if (norm <= fine_xthreshold) break;
    }
    optimize_weights_via_newton_until_converge(wstep, 0.1*fine_wthreshold, 0, max_newton_iters);
    
    m_fixed_connectivity = global_connectivity;
    return iters;
}