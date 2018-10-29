#ifndef CGAL_DEF_H
#define CGAL_DEF_H

#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Regular_triangulation_vertex_base_2.h>
#include <CGAL/Regular_triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Regular_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_2 Point_2;
typedef K::Vector_2 Vector_2;
typedef K::Segment_2 Segment_2;

// cropped voronoi
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef K::Ray_2 Ray_2;
typedef K::Line_2 Line_2;
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

typedef CGAL::Triangulation_vertex_base_with_info_2<std::pair<Point_2, double>, K> VbI;
typedef CGAL::Triangulation_face_base_with_info_2<std::vector<Point_2>, K> FbI;

typedef CGAL::Regular_triangulation_vertex_base_2<K, VbI> RT_Vb;
typedef CGAL::Regular_triangulation_face_base_2<K, FbI> RT_FbI;
typedef CGAL::Triangulation_data_structure_2<RT_Vb, RT_FbI> Tds;

typedef CGAL::Regular_triangulation_2<K, Tds> Regular_triangulation;
typedef Regular_triangulation::Bare_point BPoint;
typedef Regular_triangulation::Weighted_point WPoint;


#endif


