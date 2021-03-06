// this is adapted from src code of de Goes et al. (2012)
// Blue Noise through Optimal Transport
// from http://fernandodegoes.org/

#ifndef _PW_LINE_SEARCH_H_
#define _PW_LINE_SEARCH_H_

#include "line_search.h"

//------------//
// CWCWeights //
//------------//

template <class Scene, class T>
class CLSWeights : public CLineSearch<T, T>
{
protected:
    Scene* m_scene;
    unsigned m_nb;
    
public:
    CLSWeights(Scene* scene,
              const unsigned max_iters,
              const double max_alpha)
    : CLineSearch<T, T>(max_iters, max_alpha)
    {
        m_scene = scene;
        m_nb = m_scene->count_visible_sites();
    }
    
    double compute_function() const
    {
        return m_scene->compute_wcvt_energy();
    }
    
    void compute_gradient(std::vector<T>& V) const
    {
        m_scene->compute_weight_gradient(V);
    }
    
    bool update_scene(const std::vector<T>& X)
    {
        if (m_scene->connectivity_fixed()) return true;
        m_scene->update_weights(X, false);
        m_scene->update_triangulation();
        return has_same_vertices();
    }
    
    bool has_same_vertices() const
    {
        unsigned nb = m_scene->count_visible_sites();
        if (nb != m_nb) std::cout << "HiddenVertices: " << m_nb << " -> " << nb << std::endl;
        return (nb == m_nb);
    }
};

//-------------//
// CWCPosition //
//-------------//

template <class Scene, class Position, class Velocity>
class CLSPositions : public CLineSearch<Position, Velocity>
{
protected:
    Scene* m_scene;
    unsigned m_nb;

public:
    CLSPositions(Scene* scene,
                const unsigned max_iters,
                const double max_alpha)
    : CLineSearch<Position, Velocity>(max_iters, max_alpha)
    {
        m_scene = scene;
        m_nb = m_scene->count_visible_sites();
    }
    
    double compute_function() const
    {
        return ( - m_scene->compute_wcvt_energy() );
    }
    
    void compute_gradient(std::vector<Velocity>& V) const
    {
        m_scene->compute_position_gradient(V, -1.0);
    }
    
    bool update_scene(const std::vector<Position>& X)
    {
        if (m_scene->connectivity_fixed()) return true;
        m_scene->update_positions(X, true, false);
        m_scene->update_triangulation();
        return has_same_vertices();
    }
    
    bool has_same_vertices() const
    {
        unsigned nb = m_scene->count_visible_sites();
        if (nb != m_nb) std::cout << "HiddenVertices: " << m_nb << " -> " << nb << std::endl;
        return (nb == m_nb);
    }
};


template <class BgScene, class Position, class Velocity>
class C_BG_LSPositions : public CLineSearch<Position, Velocity>
{
protected:
	BgScene* m_scene;
	unsigned m_nb;

public:
	C_BG_LSPositions(BgScene* scene,
		const unsigned max_iters,
		const double max_alpha)
		: CLineSearch<Position, Velocity>(max_iters, max_alpha)
	{
		m_scene = scene;
		m_nb = m_scene->count_visible_sites();
	}

	double compute_function() const
	{
		return  (-m_scene->bg_compute_energy(m_scene->get_solvertype()));
	}

	void compute_gradient(std::vector<Velocity>& V) const
	{
		// -1.0 makes the result to be the gradient;
		// 1.0 results in negative gradient direction
		m_scene->bg_compute_position_gradient(V, -1.0, m_scene->get_solvertype());
	}

	bool update_scene(const std::vector<Position>& X)
	{
		if (m_scene->connectivity_fixed()) return true;
		m_scene->update_positions(X, true, false);
		m_scene->update_triangulation();
		return has_same_vertices();
	}

	bool has_same_vertices() const
	{
		unsigned nb = m_scene->count_visible_sites();
		if (nb != m_nb) std::cout << "HiddenVertices: " << m_nb << " -> " << nb << std::endl;
		return (nb == m_nb);
	}
};

#endif
