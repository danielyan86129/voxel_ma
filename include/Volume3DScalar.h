#ifndef MYVOLUME_H
#define MYVOLUME_H

#include "commondefs.h"
#include <XForm.h> // trimesh's transform matrix

//////////////////////////////////////////////////////////////////////////
/// This is the abstract class for a 3d volume storing a scalar value at each
/// voxel
//////////////////////////////////////////////////////////////////////////
class Volume3DScalar
{
public:
    Volume3DScalar() { m_sizes = ivec3(0); }
    ~Volume3DScalar() {}

    // return the value at the requested voxel
    virtual double getDataAt(int _x, int _y, int _z) const = 0;

    // return size along x-, y-, or z- dimension
    int getSizeX() const { return m_sizes[0]; }
    int getSizeY() const { return m_sizes[1]; }
    int getSizeZ() const { return m_sizes[2]; }
    const trimesh::xform& getVoxToModelMat() const { return m_vox_to_mesh; }
    const trimesh::xform& getModeltoVoxMat() const { return m_mesh_to_vox; }

protected:
    // number of voxels along x, y, and z
    ivec3 m_sizes;
    // space conversion: voxel space -> input mesh space
    trimesh::xform m_vox_to_mesh;
    // space conversion: input mesh space -> voxel space
    trimesh::xform m_mesh_to_vox;
};

#endif // MYVOLUME_H