#ifndef C_MESH_SYNTHESIZER
#define C_MESH_SYNTHESIZER

#include "algprereq.h"

#include "../DataColle/types.h"
#include "bbnot/scene.h"

#include <vector>


struct CLocalPatch
{
	RT::Vertex_handle center_vh;
	std::vector<Vector_2> nn_vecs;
	bool is_example;

	void nn_shift(int sid)
	{
		std::vector<Vector_2> cir_nn_pts = nn_vecs;
		cir_nn_pts.insert(cir_nn_pts.end(), nn_vecs.begin(), nn_vecs.end());
		int end = nn_vecs.size();
		nn_vecs.clear();
		for (int i = sid; end = sid + end; ++i)
		{
			nn_vecs.push_back(cir_nn_pts[i]);
		}
	};
};

class ALGCOLLE_CLASS CMeshSynthesizer
{
public:
	CMeshSynthesizer() {};
	~CMeshSynthesizer() {};

public:
	void set_exemplar(RT * example) { has_example_ = true; example_ = example; };
	void set_domain(Domain in_domain) { has_domain_ = true; scene_.set_domain(in_domain); };
	bool synthesize();

private:
	bool nnf_matching();

	// compute MRF patches
	void compute_local_patches(Regular_triangulation *t, std::vector<CLocalPatch> &lp_list, bool is_example);
	void get_candidate_patch(CLocalPatch &bnn, int patch_id, int start_id = 0);

	// Markovian Random Fields
	void mrf_based_capacity_assignment(CLocalPatch &slp, CLocalPatch &best_nn);
	double nnf_search(const CLocalPatch &slp, CLocalPatch &elp, int &start_id);
	double compute_matching_error(std::vector<Vector_2> s_nn_pts, std::vector<Vector_2> e_nn_pts);
	// need no shifting
	double compute_svd_matching_error(std::vector<Vector_2> s_nn_pts, std::vector<Vector_2> e_nn_pts);
	// need shifting
	double compute_index_matching_error(std::vector<Vector_2> s_nn_pts, std::vector<Vector_2> e_nn_pts);

private:
	Scene scene_;

	RT * synthss_;
	RT * example_;

	bool has_domain_ = false;
	bool has_example_ = false;

	std::vector<CLocalPatch> candidate_patches_;
	std::vector<CLocalPatch> query_patches_;
};

#endif // !C_MESH_SYNTHESIZER

