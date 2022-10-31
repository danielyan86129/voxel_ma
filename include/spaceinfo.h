#ifndef SPACE_INFO_H
#define SPACE_INFO_H

#include <cmath>
#include <memory>

#include <Cube.h> // Tao's occupancy mask
//#include <volume.h> // Tao's volume struct
#include "Volume3DScalar.h"

#include "commondefs.h"

using std::shared_ptr;

// struct SpaceInfo3D
//{
//	point
//};

struct SpaceConverter
{
    // voxel corner offset w.r.t. the center
    static const point vox_corner_offset[8];
    // return the voxel space coordinate of the requested corner
    static inline point getVoxelCorner(const ivec3& _vox, int _ix, int _iy,
                                       int _iz)
    {
        auto cidx = (_ix << 2) | (_iy << 1) | _iz;
        const auto& o = vox_corner_offset[cidx];
        return point(_vox[0] + o[0], _vox[1] + o[1], _vox[2] + o[2]);
    }
    // return whether specified voxel in the given volume is tagged as "inside"
    static bool voxTaggedAsInside(const ivec3& _vox,
                                  const shared_ptr<Volume3DScalar> _vol)
    {
        return (!SpaceConverter::outOfVolBounds(_vox, _vol)) &&
               (get_occupancy_at_vox(_vox, _vol) == 1);
    }
    // conversion from/to voxel space (continuous) to/from 3d model space
    // (continuous) (where both voxel and 3d geometry points are floating
    // numbers)
    static inline void fromVoxToModel(const point& _v, point& _p,
                                      const xform& _vox2model)
    {
        //_p = _v + point( 0.5f );
        _p = _v;
        _p = _vox2model * _p;
        // std::swap( _p[ 1 ], _p[ 0 ] );
    }
    static inline void fromModelToVox(const point& _p, point& _q,
                                      const xform& _model2vox)
    {
        // convert to voxel space
        _q = _model2vox * _p;
        // std::swap( _q[ 1 ], _q[ 0 ] );
        //_q = _q - point( 0.5f );
    }

    // conversion from/to voxel space (continuous) to/from 3d model space
    // (continuous) (where the source point is floating numbers, but target is
    // always snapped to nearest integers)
    static inline void fromVoxToModel(const point& _v, ivec3& _p,
                                      const xform& _vox2model)
    {
        // auto p = _v + point( 0.5f );
        point p;
        p = _vox2model * _v;
        _p[0] = (int)p[0];
        _p[1] = (int)p[1];
        _p[2] = (int)p[2];
        // std::swap( _p[ 0 ], _p[ 1 ] );
    }
    static inline void fromModelToVox(const point& _p, ivec3& _v,
                                      const xform& _model2vox)
    {
        // convert to voxel space
        auto q = _model2vox * _p;
        // q = q - point( 0.5f );
        // snap to nearest integer (i.e. the voxel coord.s)
        _v[0] = int(std::round(q[0]));
        _v[1] = int(std::round(q[1]));
        _v[2] = int(std::round(q[2]));
        // finally don't forget to swap x, y
        // std::swap( _v[ 1 ], _v[ 0 ] );
    }

    // test whether the given voxel is within volume bounds
    static inline bool outOfVolBounds(const ivec3& _v,
                                      const Volume3DScalar* _vol)
    {
        return (_v[0] >= _vol->getSizeX() || _v[0] < 0 ||
                _v[1] >= _vol->getSizeY() || _v[1] < 0 ||
                _v[2] >= _vol->getSizeZ() || _v[2] < 0);
    }
    static inline bool outOfVolBounds(const ivec3& _v,
                                      const shared_ptr<Volume3DScalar>& _vol)
    {
        return outOfVolBounds(_v, _vol.get());
    }

    // return the value at the given voxel. return 0 if out of bound
    static inline int get_occupancy_at_vox(const ivec3& _v,
                                           const Volume3DScalar* _vol)
    {
        auto val = SpaceConverter::outOfVolBounds(_v, _vol)
                       ? 0
                       : _vol->getDataAt(_v[0], _v[1], _v[2]) > 0.0 ? 1 : 0;
        return val;
    }
    static inline int
    get_occupancy_at_vox(const ivec3& _v,
                         const shared_ptr<Volume3DScalar>& _vol)
    {
        return get_occupancy_at_vox(_v, _vol.get());
    }

    // return the list of 8 nb voxels' coords, given a point in the voxel space
    // Note: the returned voxels are computed using constants defined in Tao's
    // geomcommon.h file
    static inline void getNbVoxelCoords(const point& _p, ivec3* _nb_vx_coords)
    {
        point v;
        // locate the given point into a voxel with lowest possible coords
        point p((int)_p[0], (int)_p[1], (int)_p[2]);
        // the eight neighbor voxels will be the offset voxels based on that one
        for (int i = 0; i < 8; ++i)
        {
            v = p + point(coordmap[i]);
            _nb_vx_coords[i][0] = (int)v[0];
            _nb_vx_coords[i][1] = (int)v[1];
            _nb_vx_coords[i][2] = (int)v[2];
        }
    }

    // return whether the nb 8 voxels of the given point are all inside the
    // volume
    static inline bool allNbVoxelsInside(const point& _p,
                                         const shared_ptr<Volume3DScalar>& _vol)
    {
        ivec3 v;
        // locate the given point into a voxel with lowest possible coords
        point p((int)_p[0], (int)_p[1], (int)_p[2]);
        // the eight neighbor voxels will be the offset voxels based on that one
        int inside_cnt = 0;
        for (int i = 0; i < 8; ++i)
        {
            v[0] = p[0] + coordmap[i][0];
            v[1] = p[1] + coordmap[i][1];
            v[2] = p[2] + coordmap[i][2];
            inside_cnt += (int)get_occupancy_at_vox(v, _vol);
        }

        return inside_cnt == 8;
    }

private:
    // constants
};

#endif // !SPACE_CONVERT_H