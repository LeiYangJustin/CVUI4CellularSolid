#ifndef C_CVT_WRAPPER_H_
#define C_CVT_WRAPPER_H_

#include "openmesh_types.h"
#include "../Optimizer/my_optimizer.h"

// generate randome double type numerics from a given interval
class WrapperUtils {
public:
	static double random_double(const double min, const double max)
	{
		double range = max - min;
		return min + (double(rand()) / double(RAND_MAX)) * range;
	}
};

class CLPCVTwrapper
{
public:
	CLPCVTwrapper() {};
	~CLPCVTwrapper() {};

	//////////
	// INIT //
	//////////

	void init(std::vector<MyMesh::Point> pts, MyMesh &mesh, unsigned int p_order = 2)
	{
		// initial
		p_order_ = p_order;
		sites_ = pts;
		set_bounding_mesh(mesh);
		build();
	}

	void init(const unsigned int nsites, const MyMesh &mesh, unsigned int p_order = 2)
	{
		// initial
		p_order_ = p_order;
		generate_random_sites(nsites, sites_);
		set_bounding_mesh(mesh);
		build();
	}

	void optimize(const unsigned int max_iter, const double end_threshold) {
		optimizer_->loop(max_iter, end_threshold);
	}

	void set_bounding_mesh(const MyMesh &in_mesh)
	{
		b_mesh_ = in_mesh; // copy mesh data
	}

private:
	unsigned int p_order_ = 2;
	double energy_;

	std::vector<MyMesh::Point> sites_;
	MyMesh b_mesh_;

	Geex::Mesh* gmesh_;

	CMyOptimizer* optimizer_;
	CLpCVT* cvt_object_;

private:
	// build
	void build()
	{
		// inti gmesh
		gmesh_ = new Geex::Mesh;
		std::vector<Geex::vec3> sites;
		for (int i = 0; i < sites_.size(); ++i)
		{
			sites.push_back(Geex::vec3(sites_[i][0], sites_[i][1], sites_[i][2]));
		}
		covert_MyMesh_to_GeexMesh(b_mesh_, *gmesh_);

		// init cvt_object
		cvt_object_ = new CLpCVT(sites, *gmesh_, p_order_);
		cvt_object_->build_with_export();

		// init optimizer
		optimizer_ = new CMyOptimizer;
		optimizer_->set_model(cvt_object_);
	};

	// generate random sites
	void generate_random_sites(const unsigned int ns, std::vector<MyMesh::Point> &pts)
	{
		double min, max;
		min = -0.2;
		max = 1.2;
		std::vector<MyMesh::Point> rand_sites;
		for (unsigned int i = 0; i < ns; ++i)
		{
			double x = WrapperUtils::random_double(min, max);
			double y = WrapperUtils::random_double(min, max);
			double z = WrapperUtils::random_double(min, max);
			rand_sites.push_back(MyMesh::Point(x, y, z));
		}
		pts = rand_sites;
	};

	// 
	void covert_MyMesh_to_GeexMesh(MyMesh mymesh, Geex::Mesh &gmesh)
	{
		std::vector<Geex::vec3> vertex_list;
		std::vector<std::vector<int>> facet_list;

		for (auto vit = mymesh.vertices_begin(); vit != mymesh.vertices_end(); ++vit)
		{
			MyMesh::Point p = mymesh.point(*vit);
			vertex_list.push_back(Geex::vec3(p[0], p[1], p[2]));
		}
		for (auto fit = mymesh.faces_begin(); fit != mymesh.faces_end(); ++fit)
		{
			std::vector<int> facet;
			for (auto fvccwit = mymesh.fv_ccwbegin(*fit); 
				fvccwit != mymesh.fv_ccwend(*fit); 
				++fvccwit)
			{
				facet.push_back(fvccwit->idx()+1);
			}
			facet_list.push_back(facet);
		}

		GeexMeshWrapper::initMesh(gmesh, vertex_list, facet_list);
	};
};

#endif // !C_LPCVT_WRAPPER_H_