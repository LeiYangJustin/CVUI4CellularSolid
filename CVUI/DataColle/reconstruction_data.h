#ifndef C_RECONSTRUCTION_DATA_H
#define C_RECONSTRUCTION_DATA_H

#include "cgal_def.h"
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/centroid.h>

class CConstraintPoint
{
public:
	CConstraintPoint() : half_width_(20), half_height_(20), point_(0.0, 0.0) {};
	CConstraintPoint(const Point_2 &_p) : half_width_(20), half_height_(20) { point_ = _p; };

	void setAllowableDirection(std::vector<Point_2> nn_pts)
	{
		// if nn_pts size == 0 -> this->point_ is a point at image border which should be fixed
		if (nn_pts.size() == 0)
			allow_dir_ = Vector_2(0, 0);
		else {
			Line_2 line;
			double good_fit = CGAL::linear_least_squares_fitting_2(nn_pts.begin(), nn_pts.end(), line, CGAL::Dimension_tag<0>());
			// if good_fit < 0.3 -> this->point_ is a juncture point for 2D ONLY
			if (good_fit < 0.3)
				allow_dir_ = Vector_2(0, 0);
			// else the fitted result is reliable
			else {
				allow_dir_ = line.to_vector();
				allow_dir_ /= sqrt(allow_dir_.squared_length());
			}
		}
	};
	const Vector_2& getAllowableDirection()
	{
		return allow_dir_;
	}

public:
	Point_2 point_;
	Vector_2 allow_dir_;
	int half_width_;
	int half_height_;
};

class CReconstructionPoint
{
public:
	CReconstructionPoint() : label_(-1), dist_(50.0), point_(0.0, 0.0) {};
	CReconstructionPoint(const Point_2 &_p, int label) : label_(label), dist_(50.0) { point_ = _p; };
	CReconstructionPoint(const Point_2 &_p) : label_(-1), dist_(50.0) { point_ = _p; };
	~CReconstructionPoint() {};
	int GetLabel() { return label_; };
	int GetDist() { return dist_; };
	void SetLabelDist(int label, double dist) { label_ = label; dist_ = dist; };

	double DistToSegment(const Segment_2 &e)
	{
		Vector_2 a = point_ - e.source();
		Vector_2 v = e.target() - e.source();
		Vector_2 proj_a_v = v*(a*v / (v*v));
		double t = sqrt(proj_a_v.squared_length()) / sqrt(v.squared_length());
		double s = a*v;
		if (s < 0)
			t = -t;
		if (t >= 0 && t <= 1)
		{
			Vector_2 d = a - proj_a_v;
			return sqrt(d.squared_length());
		}
		else if (t < 0) {
			return sqrt(a.squared_length());
		}
		else {
			return sqrt((point_ - e.target()).squared_length());
		}
	}

	struct point_less_compare {
		bool operator() (const CReconstructionPoint& lhs, const CReconstructionPoint& rhs) const {
			return lhs.point_ < rhs.point_;
		}
	};

public:
	Point_2 point_;

private:
	int label_;
	double dist_;
};


class CReconstructionEdge
{

private:
	int id_;
	Segment_2 v_edge_;
	Segment_2 f_edge_;
	std::vector<Point_2> pt_list_;
	double linearity_;
	double cos_angle_v_and_f_;
	bool good_fit_;

	double dist_thres_;
	double fabs_angle_thres_;
	double linearity_thres_;

public:
	CReconstructionEdge(): 
		id_(-1), good_fit_(false), dist_thres_(10),
		fabs_angle_thres_(0.9), linearity_thres_(0.9) {};
	~CReconstructionEdge() {};

	// id
	void set_id(const int id) {
		id_ = id;
	};
	int get_id() const {
		return id_;
	};

	// voronoi edge
	void set_voronoi_edge(const Segment_2 ve) {
		v_edge_ = ve;
	};
	const Segment_2& get_voronoi_edge() {
		return v_edge_;
	};

	// fitted edge
	void set_fitted_edge(const Segment_2 fe) {
		f_edge_ = fe;
	};
	const Segment_2& get_fitted_edge() {
		return f_edge_;
	};

	// point list
	void set_point_list(std::vector<CReconstructionPoint> pt_list)
	{
		pt_list_.clear();
		for (auto iter = pt_list.begin(); iter != pt_list.end(); ++iter)
			pt_list_.push_back(iter->point_);
		update_info();
	};

	void get_point_list(std::vector<CReconstructionPoint> &pt_list)
	{
		pt_list.clear();
		for (int i = 0; i < pt_list_.size(); i++) {
			pt_list.push_back(CReconstructionPoint(pt_list_[i], id_));
		}
			
	};

	void add_point_list(std::vector<CReconstructionPoint> pt_list)
	{
		for (auto iter = pt_list.begin(); iter != pt_list.end(); ++iter)
			pt_list_.push_back(iter->point_);
		update_info();
	};

	void add_point(CReconstructionPoint rc_pt)
	{
		pt_list_.push_back(rc_pt.point_);
		update_info();
	};

	// add a point
	bool can_add_point(CReconstructionPoint rc_pt)
	{
		std::vector<Point_2> tmp_pt_list = pt_list_;
		tmp_pt_list.push_back(rc_pt.point_);

		// fit line 
		Line_2 line;
		double tmp_linearity = linear_least_squares_fitting_2(tmp_pt_list.begin(), tmp_pt_list.end(), line, CGAL::Dimension_tag<0>());

		// finding the two ends
		Point_2 center = centroid(tmp_pt_list.begin(), tmp_pt_list.end());
		Vector_2 vdir = line.to_vector();
		double min_proj_para = 100, max_proj_para = -100;
		for (int i = 0; i < tmp_pt_list.size(); i++)
		{
			Vector_2 v1 = tmp_pt_list[i] - center;
			double tmp_proj_para = v1*vdir / sqrt(vdir.squared_length());
			if (min_proj_para > tmp_proj_para)
				min_proj_para = tmp_proj_para;
			if (max_proj_para < tmp_proj_para)
				max_proj_para = tmp_proj_para;
		}
		Segment_2 tmp_f_edge(center + min_proj_para*vdir, center + max_proj_para*vdir);
		
		// cos angle
		Vector_2 s1 = tmp_f_edge.to_vector();
		Vector_2 s2 = v_edge_.to_vector();
		double tmp_cos_angle_v_and_f = s1 *s2 / sqrt(s1.squared_length()*s2.squared_length());

		// dist
		double dist = rc_pt.DistToSegment(f_edge_);

		//
		if (tmp_linearity > std::max(linearity_thres_, 0.9*linearity_) 
			&& fabs(tmp_cos_angle_v_and_f) > std::max(fabs_angle_thres_, 0.9*tmp_cos_angle_v_and_f)
			&& dist < dist_thres_)
		{
			//std::cout << "this pt can be added because tmp_linearity(>0.3) = " << tmp_linearity 
			//	<< "cos_angle_v_and_f_(>0.7) = " << cos_angle_v_and_f_ << std::endl;
			return true;
		}
		else
		{
			//std::cout << "this pt cannot be added because tmp_linearity(>0.3) = " << tmp_linearity
			//	<< "cos_angle_v_and_f_(>0.7) = " << cos_angle_v_and_f_ << std::endl;
			return false;
		}
			
	};

	void compute_fitted_edge()
	{
		// fit line 
		Line_2 line;
		linearity_ = linear_least_squares_fitting_2(pt_list_.begin(), pt_list_.end(), line, CGAL::Dimension_tag<0>());
		// finding the two ends
		// Projection
		Point_2 center = centroid(pt_list_.begin(), pt_list_.end());
		Vector_2 vdir = line.to_vector();
		double min_proj_para = 100, max_proj_para = -100;
		for (int i = 0; i < pt_list_.size(); i++)
		{
			Vector_2 v1 = pt_list_[i] - center;
			double tmp_proj_para = v1*vdir / sqrt(vdir.squared_length());
			if (min_proj_para > tmp_proj_para)
				min_proj_para = tmp_proj_para;
			if (max_proj_para < tmp_proj_para)
				max_proj_para = tmp_proj_para;
		}
		f_edge_ = Segment_2(center + min_proj_para*vdir, center + max_proj_para*vdir);
	};

	// fitness
	double get_linearity()
	{
		// fit line 
		Line_2 line;
		linearity_ = linear_least_squares_fitting_2(pt_list_.begin(), pt_list_.end(), line, CGAL::Dimension_tag<0>());
		return linearity_;
	};

	// angle between fitted edge and voronoi edge
	double get_angle()
	{
		Vector_2 s1 = f_edge_.to_vector();
		Vector_2 s2 = v_edge_.to_vector();
		cos_angle_v_and_f_ = s1 *s2 / sqrt(s1.squared_length()*s2.squared_length());
		return cos_angle_v_and_f_;
	};

	void check_fitness() {
		//std::cout << "linearity = " << linearity_ << "; cos_angle_v_and_f_ = " << cos_angle_v_and_f_ << std::endl;
		if (linearity_ > linearity_thres_ && fabs(cos_angle_v_and_f_) > fabs_angle_thres_)
		{
			good_fit_ = true;
		}
	};

	bool is_good_fit() {
		check_fitness();
		return good_fit_;
	};

	void update_info()
	{
		compute_fitted_edge();
		get_angle();
		get_linearity();
		check_fitness();
	};
};

#endif // !C_RECONSTRUCTION_DATA_H

