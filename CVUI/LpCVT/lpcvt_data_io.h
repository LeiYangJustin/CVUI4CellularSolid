#ifndef _LPCVT_DATA_IO_
#define _LPCVT_DATA_IO_

///////////
// I / O //
///////////

#include <ostream>
#include <vector>
#include <fstream>
#include <string>

#include <LpCVT/combinatorics/delaunay.h>
#include <LpCVT/combinatorics/RVD.h>
#include <LpCVT/combinatorics/clipped_VD.h>
#include <LpCVT/common/line_stream.h>

#include "lpcvt_prereq.h"

namespace Geex {

	//==================================================================================

	/**
	* Used by save_RDT().
	*/
	class SavePrimalTriangle {
	public:
		SavePrimalTriangle(
			std::ofstream& out
		) : out_(&out) {
		}
		void operator()(unsigned int i, unsigned int j, unsigned int k) const {
			(*out_) << "f " << i + 1 << " " << j + 1 << " " << k + 1 << std::endl;
		}

	private:
		std::ofstream* out_;
	};

	/**
	* used by save_RVD()
	*/
	class SaveRVDFacets {
	public:
		SaveRVDFacets(
			std::ostream& out
		) : out_(out), cur_v_(1), cur_f_(1) {
			out << "# attribute chart facet integer" << std::endl;
		}
		void operator()(unsigned int iv, Mesh* M) const {
			for (unsigned int f = 0; f<M->nb_facets(); f++) {
				for (unsigned int i = M->facet_begin(f); i<M->facet_end(f); i++) {
					const vec3& v = M->vertex(i);
					out_ << "v " << v << std::endl;
				}
				out_ << "f ";
				for (unsigned int i = M->facet_begin(f); i<M->facet_end(f); i++) {
					out_ << cur_v_ << " ";
					cur_v_++;
				}
				out_ << std::endl;
				out_ << "# attrs f " << cur_f_ << " " << iv << std::endl;
				cur_f_++;
			}
		}
	private:
		//mutable std::ostream& out_;
		std::ostream& out_;
		mutable unsigned int cur_v_;
		mutable unsigned int cur_f_;
	};

	/**
	* Used by save_clippedVD()
	*/
	class SaveClippedVDFacets {
	public:
		SaveClippedVDFacets(
			Delaunay* delaunay, std::ostream& out, double shrink
		) : delaunay_(delaunay), out_(out), shrink_(shrink), cur_(1) {
			out << "# attribute chart facet integer" << std::endl;
		}
		void operator()(unsigned int i, int j, const vec3& p1, const vec3& p2, const vec3& p3) const {
			vec3 x0 = delaunay_->vertex(i);
			out_ << "v " << x0 + shrink_ * (p1 - x0) << std::endl;
			out_ << "v " << x0 + shrink_ * (p2 - x0) << std::endl;
			out_ << "v " << x0 + shrink_ * (p3 - x0) << std::endl;
			out_ << "f " << cur_ << " " << cur_ + 1 << " " << cur_ + 2 << std::endl;
			cur_ += 3;
			out_ << "# attrs f " << ((cur_ - 1) / 3) << " " << i << std::endl;
		}
	private:
		Delaunay* delaunay_;
		//mutable std::ostream& out_;
		std::ostream& out_;
		double shrink_;
		mutable unsigned int cur_;
	};

	//==================================================================================
}

class LPCVT_CLASS CLPCVT_IO
{
public:
	// Loads points from file
	static void load_pts(const std::string& filename, std::vector<Geex::vec3>& pts)
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

	// Loads mesh from obj file
	static void load_mesh(const std::string& filename, Geex::Mesh &m)
	{
		m.load(filename);
	}

	/**
	* Given a Restricted Voronoi Diagram, saves the Restricted Delaunay
	* Triangulation to a file in alias|wavefront .obj format.
	*/
	static void save_RDT(Geex::RestrictedVoronoiDiagram& RVD, const std::string& filename)
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

	/**
	* Saves a Restricted Voronoi Diagram to a file in alias|wavefront .obj format
	* (with Graphite extensions: facet attributes, rename as .eobj to display).
	*/
	static void save_RVD(Geex::RestrictedVoronoiDiagram& RVD, const std::string& filename)
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

	/**
	* Saves a Clipped Voronoi Diagram to a file in alias|wavefront .obj format
	* (with Graphite extensions: facet attributes, rename as .eobj to display).
	* The 'shrink' parameter is to generate the clipped Voronoi cells as in
	* Figures 4 and 5 in the paper.
	*/
	static void save_clippedVD(Geex::ClippedVoronoiDiagram & CVD, const std::string & filename, double shrink)
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

};
#endif // !_LPCVT_DATA_IO_
