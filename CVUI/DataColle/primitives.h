// this is adapted from src code of de Goes et al. (2012)
// Blue Noise through Optimal Transport
// from http://fernandodegoes.org/

#ifndef _PRIMITIVES_H_
#define _PRIMITIVES_H_

#include <vector>

#include "convex_polygon.h"

// vertex class
template <class Kernel, class Vbb>
class My_vertex_base : public Vbb
{
public:
	typedef CConvexPolygon<Kernel> ConvexPolygon;

	typedef typename Kernel::FT FT;
	
	// do not define this as Point which is a preserved key word
	// Defining this to Point will cause error related to conversion from Weighted_point to Point
	typedef typename Kernel::Point_2  Bare_Point_2; 
	typedef typename Kernel::Vector_2 Vector_2;
	typedef typename Kernel::Weighted_point_2 Weighted_point_2;

	typedef typename Vbb::Triangulation_data_structure TDS;
	typedef typename TDS::Face_handle   Face_handle;
	typedef typename TDS::Vertex_handle Vertex_handle;

	template < typename TDS2 >
	struct Rebind_TDS {
		typedef typename Vbb::template Rebind_TDS<TDS2>::Other Vb2;
		typedef My_vertex_base<Kernel, Vb2> Other;
	};

public:
	My_vertex_base() : Vbb()
	{
		index_ = -1;
	}
	My_vertex_base(const Weighted_point_2 & p) : Vbb(p)
	{
		index_ = -1;
	}
	My_vertex_base(const Weighted_point_2 & p, Face_handle f) : Vbb(p, f)
	{
		index_ = -1;
	}

	// ID //
	int get_index() const
	{
		return this->index_;
	}
	void set_index(const int &id)
	{
		index_ = id;
	}

	// POSITION / WEIGHT / CAPACITY //
	const Bare_Point_2& get_position() const
	{
		return this->point().point();
	}
	void set_position(const Bare_Point_2& p)
	{
		FT w = get_weight();
		Weighted_point wp(p, w);
		this->set_point(wp);
	}
	const FT get_weight() const
	{
		return this->point().weight();
	}
	void set_weight(const FT w)
	{
		const Bare_Point_2& p = get_position();
		Weighted_point wp(p, w);
		this->set_point(wp);
	}
	const FT get_capacity() const
	{
		return this->capacity_;
	}
	void set_capacity(const FT cap)
	{
		this->capacity_ = cap;
	}

	// CELL //
	void reset_cell()
	{
		cell_.clear();
	}
	void set_cell(const ConvexPolygon& shape)
	{
		cell_ = shape;
	}

	// AREA ; VARIANCE ; CENTROID //
	FT compute_area() const
	{
		return cell_.compute_area();
	}
	Bare_Point_2 compute_centroid() const
	{
		return cell_.compute_centroid();
	}
	FT compute_variance() const
	{
		return cell_.compute_variance(get_position());
	}

private:
	int index_;
	ConvexPolygon cell_;
	double capacity_; // constraint
};


template <class Kernel, class Fbb>
class My_face_base : public Fbb
{
public:
	typedef typename Fbb::Triangulation_data_structure  TDS;
	typedef typename TDS::Vertex_handle Vertex_handle;
	typedef typename TDS::Face_handle   Face_handle;

	template < typename TDS2 >
	struct Rebind_TDS {
		typedef typename Fbb::template Rebind_TDS<TDS2>::Other Fb2;
		typedef My_face_base<Kernel, Fb2> Other;
	};

public:
	My_face_base()
		: Fbb()
	{
	}

	My_face_base(Vertex_handle v1,
		Vertex_handle v2,
		Vertex_handle v3)
		: Fbb(v1, v2, v3)
	{
	}

	My_face_base(Vertex_handle v1,
		Vertex_handle v2,
		Vertex_handle v3,
		Face_handle f1,
		Face_handle f2,
		Face_handle f3)
		: Fbb(v1, v2, v3, f1, f2, f3)
	{
	}
};

#endif // _PRIMITIVES_H_
