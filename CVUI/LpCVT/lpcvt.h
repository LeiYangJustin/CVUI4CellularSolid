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

#ifndef C_LPCVT_H_
#define C_LPCVT_H_

#include <ostream>
#include <vector>
#include <fstream>
#include <string>

#include <LpCVT/combinatorics/delaunay.h>
#include <LpCVT/combinatorics/RVD.h>
#include <LpCVT/combinatorics/clipped_VD.h>
#include <LpCVT/algebra/F_Lp.h>
#include <LpCVT/common/line_stream.h>

#include "lpcvt_prereq.h"
#include "lpcvt_data_io.h"

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

class LPCVT_CLASS GeexMeshWrapper 
{
public:
	static void initMesh(Geex::Mesh &gmesh,
		const std::vector<Geex::vec3> vertex_list,
		const std::vector<std::vector<int>> facet_list)
	{
		gmesh.init(vertex_list, facet_list);
	}
};

class LPCVT_CLASS CLpCVT
{
private:
	// this would be changed
	Geex::Mesh* b_mesh_;
	
	std::vector<Geex::vec3> sites_;
	std::vector<Geex::vec3> gradients_;

	double energy_;
	unsigned int p_order_;

	bool is_export_=false;

public:
	CLpCVT() {};
	CLpCVT(std::vector<Geex::vec3> pts, Geex::Mesh &in_mesh, int p_order) {
		p_order_ = p_order;
		sites_ = pts;
		b_mesh_ = &in_mesh;
	};
	~CLpCVT() {};

	///////////////
	// GET / SET //
	///////////////

	void get_sites(std::vector<Geex::vec3> &pts)
	{
		pts = sites_;
	};

	void set_sites(const std::vector<Geex::vec3> &pts)
	{
		sites_ = pts;
	}

	////////////////
	// Attributes //
	////////////////

	void cvt_gradient(std::vector<Geex::vec3> &gradients) {
		gradients = gradients_;
	};

	double cvt_energy() {
		return energy_;
	};

	////////////
	// Update //
	////////////

	void update_decomposition(const std::vector<Geex::vec3> &pts, bool is_export = false)
	{
		clean_cvt();
		set_sites(pts);

		if (!is_export)
			build();
		else
			build_with_export();
	};

	// build the decomposition and compute the objective function and gradients
	void build()
	{
		is_export_ = false;
		assert(sites_.size() != 0);
		compute_F_g(energy_, gradients_);
	};

	void build_with_export()
	{
		is_export_ = true;
		assert(sites_.size() != 0);
		compute_F_g(energy_, gradients_);
	};

private:
	void get_combinatorics(
		std::vector<int>& I, std::vector<Geex::vec3>& C, std::vector<int>& F, bool volume
	);

	void clean_cvt()
	{
		sites_.clear();
		gradients_.clear();
	}

	// compute objective function and ist gradient vector
	void compute_F_g();
	void compute_F_g(double &energy, std::vector<Geex::vec3> &gradients);
};

#endif // !C_LPCVT_H_

