// this is adapted from src code of de Goes et al. (2012)
// Blue Noise through Optimal Transport
// from http://fernandodegoes.org/

#ifndef _SCENE_H_
#define _SCENE_H_

// STL
#include <vector>

// EIGEN
#include <Eigen/Dense>
#include <Eigen/Sparse>

// local
#include "../../DataColle/types.h"
#include "../../DataColle/util.h"

class Scene
{    
private:
    RT m_rt;
    Domain m_domain;
    std::vector<FT> m_capacities;
    std::vector<Vertex_handle> m_vertices;
    bool m_fixed_connectivity;
    double m_tau;
    
    //bool m_timer_on;
    //std::vector<double> m_timer;

	double M_PI = 3.14159265359;

public:
    Scene()
    {
        srand(0);
        m_tau = 1.0;
        init_domain(4);
        m_rt.set_domain(&m_domain);
        m_fixed_connectivity = false;
        //m_timer_on = false;
    }
    
    ~Scene()
    {
        clear();
    }    
    
	// ADDED INTERFACES //
	// interface to set a user-prescribed domain
	void set_domain(Domain in_domain) { 
		m_domain.clear();
		m_domain = in_domain; 
		m_rt.set_domain(&m_domain); 
	};

	// interface between matching and optimization
	void set_rt(RT * rt) { m_rt = *rt; };
	void get_rt(RT * rt) { rt = &m_rt; };

	// optimize with specific constraint
	unsigned optimize_all_with_background_constraint(FT wstep, FT xstep, 
		unsigned max_newton_iters,
		FT epsilon, unsigned max_iters,
		std::ostream& out);

	// ORIGINAL CODES //

    double get_tau() const { return m_tau; }
    void set_tau(double tau) { m_tau = tau; }

    //void toggle_timer() { m_timer_on = !m_timer_on; }
    
    void toggle_connectivity() { m_fixed_connectivity = !m_fixed_connectivity; }
    
    bool connectivity_fixed() const { return m_fixed_connectivity; }
    
    void clear()
    {
        clear_triangulation();
    }

    // IO //
    
    void init_domain(unsigned nb);
        
    //unsigned load(const QString& filename);
    //
    //void load_dat_file(const QString& filename, std::vector<Point_2>& points) const;
    //
    //void load_txt_file(const QString& filename, std::vector<Point_2>& points) const;

    //void save(const QString& filename) const;
    //
    //void save_dat(const QString& filename, const std::vector<Point_2>& points) const;
    //
    //void save_txt(const QString& filename, const std::vector<Point_2>& points) const;
    //
    //void save_eps(const QString& filename) const;
    
	bool Scene::write_updated_triangulation(std::string fname, RT &rt);

    // SITES //
    
    void generate_random_sites(const unsigned nb);

    void generate_regular_grid(const unsigned nx, const unsigned ny);
    
    void generate_hextille_grid(const unsigned nb);
    
    // RENDER //
    
    //void draw_point(const Point_2& a) const; 
    //
    //void draw_segment(const Point_2& a, const Point_2& b) const;    
    //
    //void draw_triangle(const Point_2& a, const Point_2& b, const Point_2& c) const;    
    //
    //void draw_polygon(const std::vector<Point_2>& polygon) const;
    //
    //void draw_region(const std::vector<Point_2>& polygon) const;
    //
    //void draw_circle(const Point_2& center, 
    //                 const FT scale, 
    //                 const std::vector<Point_2>& pts) const;

    //void draw_domain(const float line_width,
    //                 const float red,
    //                 const float green,
    //                 const float blue) const;
    //
    //void draw_sites(const float point_size,
    //                const float red,
    //                const float green,
    //                const float blue) const;
    //
    //void draw_vertices(const float point_size) const;
    //
    //void draw_faces(const float red, 
    //                const float green,
    //                const float blue) const;
    //
    //void draw_primal(const float line_width,
    //                 const float red,
    //                 const float green,
    //                 const float blue) const;
    //
    //void draw_dual(const float line_width,
    //               const float red,
    //               const float green,
    //               const float blue) const;

    //void draw_bounded_dual(const float line_width,
    //                       const float red,
    //                       const float green,
    //                       const float blue) const;
    //
    //void draw_weights() const;
    //
    //void draw_capacity() const;
    //
    //void draw_variance() const;

    //void draw_regular_sites() const;

    //void draw_regularity() const;

    //void draw_cell(Vertex_handle vertex, bool filled) const;
    //
    //void draw_capacity_histogram(const unsigned nbins, 
    //                             const double xmin,
    //                             const double xmax,
    //                             const double ymin,
    //                             const double ymax) const;
    //
    //void draw_weight_histogram(const double range,
    //                           const unsigned nbins, 
    //                           const double xmin,
    //                           const double xmax,
    //                           const double ymin,
    //                           const double ymax) const;
    //
    //void draw_histogram(const std::vector<unsigned>& histogram,
    //                    const double xmin,
    //                    const double xmax,
    //                    const double ymin,
    //                    const double ymax) const;
    
    // HISTOGRAM //
    
    void compute_capacity_histogram(std::vector<unsigned>& histogram) const;    
    
    void compute_weight_histogram(const double range, std::vector<unsigned>& histogram) const;    
    
    // INIT //
    
    FT compute_mean_capacity() const;

    unsigned count_visible_sites() const;
    
    void collect_visible_points(std::vector<Point_2>& points) const;
    
    void collect_visible_weights(std::vector<FT>& weights) const;
    
    void collect_sites(std::vector<Point_2>& points,
                       std::vector<FT>& weights) const;

    void clear_triangulation();
    
    void update_triangulation();
    
    void construct_triangulation(const std::vector<Point_2>& points,
                                 const std::vector<FT>& weights);
    
    void populate_vertices(const std::vector<Point_2>& points,
                           const std::vector<FT>& weights);
    
    Vertex_handle insert_vertex(const Point_2& point,
                                const FT weight,
                                const unsigned index);
    
    void compute_capacities(std::vector<FT>& capacities) const;

    void reset_weights();
    
    void update_positions(const std::vector<Point_2>& points, bool clamp = true, bool hidden = true);
                          
    void update_weights(const std::vector<FT>& weights, bool hidden = true);
    
    void project_positions_to_domain();
    
    // ENERGY //
    
    FT compute_wcvt_energy();
    
    FT compute_total_area() const;
    
    void compute_position_gradient(std::vector<Vector_2>& gradient, FT coef = 1.0);
    
    void compute_weight_gradient(std::vector<FT>& gradient, FT coef = 1.0);
    
    FT compute_weight_threshold(FT epsilon) const;

    FT compute_position_threshold(FT epsilon) const;

    // OPTIMIZER //
    
    FT optimize_positions_via_lloyd(bool update);
    
    FT optimize_positions_via_gradient_ascent(FT timestep, 
                                              bool update);

    FT optimize_weights_via_gradient_descent(FT timestep, 
                                             bool update);    
    
    FT optimize_weights_via_newton(FT timestep,
                                   bool update);
    
    unsigned optimize_weights_via_gradient_descent_until_converge(FT timestep, 
                                                                  FT threshold,
                                                                  unsigned update,
                                                                  unsigned max_iters);
    
    unsigned optimize_weights_via_newton_until_converge(FT timestep, 
                                                        FT threshold,
                                                        unsigned update,
                                                        unsigned max_iters);
    
    unsigned optimize_all(FT wstep, FT xstep, unsigned max_newton_iters,
                          FT epsilon, unsigned max_iters,
                          std::ostream& out);

    bool solve_newton_step(const std::vector<FT>& b, 
                           std::vector<FT>& x);
    
    void build_laplacian(const FT scale,
                         const std::map<unsigned, unsigned>& indices,
                         Eigen::SparseMatrix<double>& A) const;
    
    bool solve_linear_system(const Eigen::SparseMatrix<double>& A,
                             std::vector<double>& x,
                             const std::vector<double>& b) const;

    // REGULARITY //
    
    FT compute_regularity_threshold() const;
    
    void detect_and_break_regularity(FT max_radius, unsigned max_teleport);

    FT compute_max_regularity() const;
    
    void compute_variance_vector(std::vector<FT>& variance) const;
    
    void compute_regularity_vector(const std::vector<FT>& variance,
                                   std::vector<FT>& regularity) const;
    
    FT compute_regularity(Vertex_handle vi, const std::vector<FT>& variance) const;
    
    void jitter_vertices(const std::set<Vertex_handle>& vertices, const FT max_radius);
    
    Point_2 jitter_point(const Point_2& p, const FT max_radius) const;
};

#endif
