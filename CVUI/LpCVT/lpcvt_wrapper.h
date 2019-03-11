// This file is adopted from LpCVT main.cpp

// Quick reference:
//
// Mesh: stores a Piecewise Linear Complex (+ functions to read it from a file in .obj format)
// Delaunay: abstract interface to Delaunay triangulation (+ Delaunay_CGAL: implementation with CGAL)
// RestrictedVoronoiDiagram: Given a Mesh and a Delaunay, computes the intersection between a Voronoi diagram and a Mesh. 
// ClippedVoronoiDiagram: computes the intersection between a Voronoi diagram and the interior of a closed Mesh
//
// RestrictedVoronoiDiagram and ClippedVoronoiDiagram are accessed through the for_each_triangle(do_it) function,
// that does the actual computation, and calls the user-provided function do_it(i,j,k,l) for each integration
// simplex. Note that do_it can be a function or an object that overloards operator(). 
// This mechanism can be used for both computing F-Lp and displaying the clipped Voronoi cells as in figures 4,5,6
// in the paper.

#ifndef C_LPCVT_WRAPPER_H_
#define C_LPCVT_WRAPPER_H_

#include "lpcvt_prereq.h"

#include <ostream>
#include <vector>
#include <fstream>
#include <string>

#include <LpCVT/combinatorics/delaunay.h>
#include <LpCVT/combinatorics/RVD.h>
#include <LpCVT/combinatorics/clipped_VD.h>
#include <LpCVT/algebra/F_Lp.h>
#include <LpCVT/common/line_stream.h>

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

namespace Geex {

	/**
	* Used by get_combinatorics() in volume mode
	*/
	class MemorizeIndices {
	public:
		MemorizeIndices(
			std::vector<int>& I_in,
			std::vector<vec3>& C_in
		) : I(I_in), C(C_in) {
			I.resize(0);
			C.resize(0);
		}

		void operator() (
			unsigned int i,
			int j,
			const VertexEdge& v1,
			const VertexEdge& v2,
			const VertexEdge& v3
			) const {
			I.push_back(i);
			I.push_back(v1.sym[2]);
			I.push_back(v1.sym[1]);
			I.push_back(v1.sym[0]);
			I.push_back(v2.sym[2]);
			I.push_back(v2.sym[1]);
			I.push_back(v2.sym[0]);
			I.push_back(v3.sym[2]);
			I.push_back(v3.sym[1]);
			I.push_back(v3.sym[0]);
			C.push_back(v1);
			C.push_back(v2);
			C.push_back(v3);
		}
	private:
		//mutable std::vector<int>& I;
		std::vector<int>& I;
		//mutable std::vector<vec3>& C;
		std::vector<vec3>& C;
	};

	/**
	* Used by get_combinatorics() in surface mode
	*/
	class MemorizeIndicesAndFacets {
	public:
		MemorizeIndicesAndFacets(
			const RestrictedVoronoiDiagram& RVD_in,
			std::vector<int>& I_in,
			std::vector<vec3>& C_in,
			std::vector<int>& F_in
		) : RVD(RVD_in), I(I_in), C(C_in), F(F_in) {
			I.resize(0);
			C.resize(0);
			F.resize(0);
		}

		void operator() (
			unsigned int i,
			const VertexEdge& v1,
			const VertexEdge& v2,
			const VertexEdge& v3
			) const {
			I.push_back(i);
			I.push_back(v1.sym[2]);
			I.push_back(v1.sym[1]);
			I.push_back(v1.sym[0]);
			I.push_back(v2.sym[2]);
			I.push_back(v2.sym[1]);
			I.push_back(v2.sym[0]);
			I.push_back(v3.sym[2]);
			I.push_back(v3.sym[1]);
			I.push_back(v3.sym[0]);
			F.push_back(RVD.current_facet());
			C.push_back(v1);
			C.push_back(v2);
			C.push_back(v3);
		}
	private:
		const RestrictedVoronoiDiagram& RVD;
		//mutable std::vector<int>& I;
		//mutable std::vector<vec3>& C;
		//mutable std::vector<int>& F;
		std::vector<int>& I;
		std::vector<vec3>& C;
		std::vector<int>& F;
	};
}

class LPCVT_CLASS CLPCVTwrapper
{
public:
	CLPCVTwrapper();
	~CLPCVTwrapper();

	void test_all()
	{
		//
		std::string mesh_filename_ = "D:\\MyProjects\\CVUI4CellularSolid\\CVUI\\MainPrj3D\\test_data\\three_holes.obj";
		std::string pts_filename_ = "D:\\MyProjects\\CVUI4CellularSolid\\CVUI\\MainPrj3D\\test_data\\three_holes.pts";
		
		//
		std::cerr << "============= geometry->combinatorics test ==========" << std::endl;
		test_combinatorics(mesh_filename_, pts_filename_);
		std::cerr << "============= combinatorics->algebra test  ==========" << std::endl;
		std::cerr << "(note: expect large values for f and g)" << std::endl;
		test_algebra(mesh_filename_, pts_filename_);
	}

	void test_combinatorics(const std::string& mesh_filename, const std::string& pts_filename);

	void test_algebra(const std::string& mesh_filename, const std::string& pts_filename);

private:

	/**
	* Loads points from file
	*/
	void load_pts(const std::string& filename, std::vector<Geex::vec3>& pts);

	/**
	* Given a Restricted Voronoi Diagram, saves the Restricted Delaunay
	* Triangulation to a file in alias|wavefront .obj format.
	*/
	void save_RDT(Geex::RestrictedVoronoiDiagram& RVD, const std::string& filename);

	/**
	* Saves a Restricted Voronoi Diagram to a file in alias|wavefront .obj format
	* (with Graphite extensions: facet attributes, rename as .eobj to display).
	*/
	void save_RVD(Geex::RestrictedVoronoiDiagram& RVD, const std::string& filename);


	/**
	* Saves a Clipped Voronoi Diagram to a file in alias|wavefront .obj format
	* (with Graphite extensions: facet attributes, rename as .eobj to display).
	* The 'shrink' parameter is to generate the clipped Voronoi cells as in
	* Figures 4 and 5 in the paper.
	*/
	void save_clippedVD(Geex::ClippedVoronoiDiagram& CVD, const std::string& filename, double shrink);
	
	//===================================================================================================
	//
	// Geometry and combinatorics tests
	//
	//  Computes and output to a file:
	//     Restricted Voronoi Diagram (intersection between a mesh and a 3D Voronoi diagram)
	//     Restricted Delaunay Triangulation (dual of a Restricted Voronoi Diagram)
	//     Clipped Voronoi Diagram (intersection between the interior of a mesh and a 3D Voronoi diagram)
	//
	//=================================================================================================== 

	// ALGEBRA
	//===================================================================================================
	//
	// Algebra tests
	//
	//  Computes :
	//     F_{L_p} in the surface case
	//     F_{L_p} in the volume case
	//
	//=================================================================================================== 

	/**
	* Gets the combinatorics of the integration simplices,
	* i.e. 10 integers per integration simplex.
	* (see Section 3.1 in the paper)
	* Returns also the array of C vertices (three per integration simplex).
	* Since they are easy to get during the combinatorial phase, they are
	* computed here and kept for the algebraic phase.
	*
	* In 2D mode (volume = false), returns also the array F.
	*   F[i] indicates the facet that contains the i-th integration simplex.
	*
	*/
	void get_combinatorics(
		Geex::Mesh* M, const std::vector<Geex::vec3>& pts,
		std::vector<int>& I, std::vector<Geex::vec3>& C, std::vector<int>& F, bool volume
	); 
	
	/**
	* Computes F_{L_p} and its gradient.
	*/
	void compute_F_g(Geex::Mesh* m, const std::vector<Geex::vec3>& pts, unsigned int p, bool volume);
};


#endif // !C_LPCVT_WRAPPER_H_