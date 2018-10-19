#ifndef CGAL_DEF_H
#define CGAL_DEF_H

#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Regular_triangulation_vertex_base_2.h>
#include <CGAL/Regular_triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Regular_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Point_2<K> Point2;
typedef CGAL::Vector_2<K> Vector2;
typedef CGAL::Segment_2<K> Segment2;

//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//
//typedef My_vertex_base<K> Vb;
//typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
//typedef CGAL::Triangulation_2<K, Tds> Triangulation;

//struct LabelDist
//{
//	int label_;
//	double dist_;
//};

typedef CGAL::Triangulation_vertex_base_with_info_2<std::vector<Point2>, K> VbI;
typedef CGAL::Regular_triangulation_vertex_base_2<K, VbI> RT_VbI;
typedef CGAL::Regular_triangulation_face_base_2<K> RT_Fb;
typedef CGAL::Triangulation_data_structure_2<RT_VbI, RT_Fb> Tds;

typedef CGAL::Regular_triangulation_2<K, Tds> Regular_triangulation;
typedef Regular_triangulation::Bare_point BPoint;
typedef Regular_triangulation::Weighted_point WPoint;

#endif


