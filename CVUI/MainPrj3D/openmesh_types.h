#ifndef C_MY_OPENMESH_TYPES_H_
#define C_MY_OPENMESH_TYPES_H_

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

struct MyTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::Vec3d Point; // use double-values points
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  MyMesh;

#endif // !C_MY_OPENMESH_TYPES_H_
