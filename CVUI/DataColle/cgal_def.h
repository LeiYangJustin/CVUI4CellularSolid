#ifndef CGAL_DEF_H
#define CGAL_DEF_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_2<K> Regular_triangulation;
typedef Regular_triangulation::Bare_point BPoint;
typedef Regular_triangulation::Weighted_point WPoint;

#endif


