#include "lpcvt.h"

void CLpCVT::get_combinatorics(std::vector<int>& I, std::vector<Geex::vec3>& C, std::vector<int>& F, bool volume)
{
	//
	using namespace Geex;

	//
	Delaunay* delaunay = Delaunay::create("CGAL");
	delaunay->set_vertices(sites_);
	if (volume) {
		ClippedVoronoiDiagram CVD(delaunay, b_mesh_);
		CVD.for_each_triangle(MemorizeIndices(I, C));

		if (is_export_) CLPCVT_IO::save_clippedVD(CVD, "cvd.obj", 0.7);
	}
	else {
		RestrictedVoronoiDiagram RVD(delaunay, b_mesh_);
		RVD.set_symbolic(true);
		RVD.for_each_triangle(MemorizeIndicesAndFacets(RVD, I, C, F));
		
		if (is_export_) {
			CLPCVT_IO::save_RVD(RVD, "rvd.obj");
			CLPCVT_IO::save_RDT(RVD, "rdt.obj");
		}
	}

	delete delaunay;
}

void CLpCVT::compute_F_g()
{
	//
	using namespace Geex;

	//
	std::cerr << "nb pts = " << sites_.size() << "   nb facets = " << b_mesh_->nb_facets() << std::endl;
	std::vector<int> I;
	std::vector<vec3> C;
	std::vector<int> F;
	bool volume = true;
	get_combinatorics(I, C, F, volume);
	unsigned int nb_integration_simplices = (unsigned int)I.size() / 10;
	std::vector<mat3> M(nb_integration_simplices);
	for (unsigned int i = 0; i<M.size(); i++) {
		M[i].load_identity();
		// or replace with anisotropy field
		//   In 2D: use F[i] to retreive the index of the facet that contains
		//      the current integration simplex (and access an array of per-facet anisotropy).
		//   In 3D: use geometric search from the centroid of the current
		//      integration simplex.
	}
	std::vector<plane3> Q(b_mesh_->nb_facets());
	for (unsigned int i = 0; i<b_mesh_->nb_facets(); i++) {
		Q[i] = b_mesh_->facet_plane(i);
	}
	std::vector<double> g(sites_.size() * 3);
	double f = compute_F_Lp(volume, p_order_, b_mesh_, I, C, sites_, Q, M, g);

	// output
	energy_ = f;
	gradients_.clear();
	for (unsigned int i = 0; i < sites_.size(); ++i) {
		vec3 v(g[i * 3], g[i * 3 + 1], g[i * 3 + 2]);
		gradients_.push_back(v);
	}
}

void CLpCVT::compute_F_g(double &energy, std::vector<Geex::vec3> &gradients)
{
	compute_F_g();
	energy = energy_;
	gradients = gradients_;
}