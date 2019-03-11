#include "lpcvt_wrapper.h"

CLPCVTwrapper::CLPCVTwrapper()
{
}


CLPCVTwrapper::~CLPCVTwrapper()
{
}

void CLPCVTwrapper::test_combinatorics(
	const std::string & mesh_filename, 
	const std::string & pts_filename)
{
	//
	using namespace Geex;

	//
	Mesh M;
	unsigned int nb_borders = M.load(mesh_filename);
	std::vector<vec3> pts;
	load_pts(pts_filename, pts);
	Delaunay* delaunay = Delaunay::create("CGAL");
	RestrictedVoronoiDiagram RVD(delaunay, &M);
	ClippedVoronoiDiagram CVD(delaunay, &M);

	delaunay->set_vertices(pts);
	save_RVD(RVD, "rvd.obj");
	save_RDT(RVD, "rdt.obj");
	if (nb_borders == 0) {
		save_clippedVD(CVD, "cvd.obj", 0.7);
	}
	delete delaunay;
}

void CLPCVTwrapper::test_algebra(
	const std::string & mesh_filename, 
	const std::string & pts_filename)
{
	//
	using namespace Geex;

	//
	Mesh M;
	unsigned int nb_borders = M.load(mesh_filename);
	std::vector<vec3> pts;
	load_pts(pts_filename, pts);
	if (nb_borders == 0) {
		std::cerr << "          ========== volume LpCVT test ======" << std::endl;
		compute_F_g(&M, pts, 4, true);
	}
	std::cerr << "          ========== surface LpCVT test ======" << std::endl;
	compute_F_g(&M, pts, 4, false);
}

void CLPCVTwrapper::load_pts(const std::string & filename, std::vector<Geex::vec3>& pts)
{
	//
	using namespace Geex;

	//
	pts.clear();
	std::ifstream in_stream(filename.c_str());
	if (!in_stream) {
		std::cerr << "Could not open " << filename << std::endl;
		return;
	}
	LineInputStream in(in_stream);
	while (!in.eof()) {
		in.get_line();
		std::string kw;
		in >> kw;
		if (kw == "v") {
			vec3 v;
			in >> v;
			pts.push_back(v);
		}
	}
}

void CLPCVTwrapper::save_RDT(
	Geex::RestrictedVoronoiDiagram & RVD, 
	const std::string & filename)
{
	//
	using namespace Geex;

	//
	std::cerr << "Computing and saving RDT to " << filename << std::endl;
	std::ofstream out(filename.c_str());
	for (unsigned int i = 0; i<RVD.delaunay()->nb_vertices(); i++) {
		out << "v " << RVD.delaunay()->vertex(i) << std::endl;
	}
	RVD.for_each_primal_triangle(SavePrimalTriangle(out));
	out.close();
	std::cerr << "Done." << std::endl;
}

void CLPCVTwrapper::save_RVD(
	Geex::RestrictedVoronoiDiagram & RVD, 
	const std::string & filename)
{
	//
	using namespace Geex;

	//
	std::ofstream out(filename.c_str());
	if (!out) {
		std::cerr << "could not open file." << std::endl;
		return;
	}
	std::cerr << "Computing and saving RVD" << std::endl;
	bool sym_backup = RVD.symbolic();
	RVD.set_symbolic(true);
	RVD.for_each_facet(SaveRVDFacets(out));
	RVD.set_symbolic(sym_backup);
	std::cerr << "Saved RVD in " << filename << std::endl;
}

void CLPCVTwrapper::save_clippedVD(
	Geex::ClippedVoronoiDiagram & CVD, 
	const std::string & filename, 
	double shrink)
{
	//
	using namespace Geex;

	//
	std::ofstream out(filename.c_str());
	if (!out) {
		std::cerr << "could not open file." << std::endl;
		return;
	}
	std::cerr << "Computing and saving clipped VD" << std::endl;
	CVD.for_each_triangle(SaveClippedVDFacets(CVD.delaunay(), out, shrink));
	std::cerr << "Saved clipped VD in " << filename << std::endl;
}

void CLPCVTwrapper::get_combinatorics(
	Geex::Mesh * M, 
	const std::vector<Geex::vec3>& pts, 
	std::vector<int>& I, 
	std::vector<Geex::vec3>& C, 
	std::vector<int>& F, 
	bool volume)
{
	//
	using namespace Geex;

	//
	Delaunay* delaunay = Delaunay::create("CGAL");
	delaunay->set_vertices(pts);
	if (volume) {
		ClippedVoronoiDiagram CVD(delaunay, M);
		CVD.for_each_triangle(MemorizeIndices(I, C));
	}
	else {
		RestrictedVoronoiDiagram RVD(delaunay, M);
		RVD.set_symbolic(true);
		RVD.for_each_triangle(MemorizeIndicesAndFacets(RVD, I, C, F));
	}
	delete delaunay;
}

void CLPCVTwrapper::compute_F_g(Geex::Mesh * m, const std::vector<Geex::vec3>& pts, unsigned int p, bool volume)
{
	//
	using namespace Geex;

	//
	std::cerr << "nb pts = " << pts.size() << "   nb facets = " << m->nb_facets() << std::endl;
	std::vector<int> I;
	std::vector<vec3> C;
	std::vector<int> F;
	get_combinatorics(m, pts, I, C, F, volume);
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
	std::vector<plane3> Q(m->nb_facets());
	for (unsigned int i = 0; i<m->nb_facets(); i++) {
		Q[i] = m->facet_plane(i);
	}
	std::vector<double> g(pts.size() * 3);
	double f = compute_F_Lp(volume, p, m, I, C, pts, Q, M, g);
	double gnorm = 0.0;
	for (unsigned int i = 0; i<g.size(); i++) {
		gnorm += g[i] * g[i];
	}
	gnorm = ::sqrt(gnorm);
	std::cerr.precision(16);
	std::cerr << (volume ? "volume " : "surface ")
		<< "F_L" << p << ":"
		<< "f=" << std::scientific << f << "  g=" << gnorm << std::endl;
}

