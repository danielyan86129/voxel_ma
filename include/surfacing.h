#ifndef SURFACING_H
#define SURFACING_H

#include <array>
#include <memory>
#include <set>
#include <vector>

// #include <volume.h> // Tao's volume struct
#include "Volume3DScalar.h"

#include "commondefs.h"
#include "octree.h"
#include "reterrcode.h"
#include "spaceinfo.h"

enum class SurfacerErrCode
{
    // failure
    FAILURE,
    // success
    SUCCESS,
    // err: empty volume was given
    EMPTY_VOL,
    // err: volume type not supported by surfacer
    INVALID_VOL_TYPE
};

using std::array;
using std::set;
using std::shared_ptr;
using std::vector;

class Surfacer
{
public:
    /*Constructors*/
    Surfacer();
    // constructs from a volume struct
    Surfacer(shared_ptr<Volume3DScalar> _vol);
    /*Destructor*/
    ~Surfacer();

    /*major interfaces*/
    // extract boundary vertices and faces from the volume.
    // each face will be a quad (4-tuple of indices into the returned vts list)
    SurfacerErrCode
    extractBoundaryVtsAndFaces(const shared_ptr<Volume3DScalar>& _vol,
                               vector<point>& _vts, vector<ivec4>& _faces,
                               bool _need_euler);
    // each face will be a triangle
    SurfacerErrCode
    extractBoundaryVtsAndFaces(const shared_ptr<Volume3DScalar>& _vol,
                               vector<point>& _vts, vector<uTriFace>& _faces,
                               eulerchar* _euler = nullptr);
    // extract only boundary vts from the volume
    SurfacerErrCode extractBoundaryVts(const shared_ptr<Volume3DScalar>& _vol,
                                       vector<point>& _vts);
    // extract boundary for sparse-space-partition structure
    SurfacerErrCode
    extractBoundaryVtsAndFaces(const shared_ptr<OctreeVolume>& _octree,
                               vector<point>& _vts, vector<ivec4>& _faces,
                               bool _need_euler);
    // @return the euler characteristic.
    // @pre face must be 4-tuple of vert indices that trace the boundary if
    // followed in order
    int computeEulerChar(const vector<point>& _vts, const vector<ivec4>& _faces,
                         const shared_ptr<Volume3DScalar>& _vol,
                         eulerchar& _euler) const;

protected:
    /*helpers*/
    // @return 8 corners' coords at voxel i,j,k. Note: must be within volume
    inline void get_corners_coords(int _i, int _j, int _k,
                                   const shared_ptr<Volume3DScalar>& _vol,
                                   ivec3* _crns)
    {
        // const auto& mat = _vol->getVoxToModelMat();
        auto mat = trimesh::xform::identity();
        SpaceConverter::fromVoxToModel(point(_i + 0.5f, _j - 0.5f, _k - 0.5f),
                                       _crns[0], mat);
        SpaceConverter::fromVoxToModel(point(_i + 0.5f, _j + 0.5f, _k - 0.5f),
                                       _crns[1], mat);
        SpaceConverter::fromVoxToModel(point(_i - 0.5f, _j - 0.5f, _k - 0.5f),
                                       _crns[2], mat);
        SpaceConverter::fromVoxToModel(point(_i - 0.5f, _j + 0.5f, _k - 0.5f),
                                       _crns[3], mat);
        SpaceConverter::fromVoxToModel(point(_i + 0.5f, _j - 0.5f, _k + 0.5f),
                                       _crns[4], mat);
        SpaceConverter::fromVoxToModel(point(_i + 0.5f, _j + 0.5f, _k + 0.5f),
                                       _crns[5], mat);
        SpaceConverter::fromVoxToModel(point(_i - 0.5f, _j - 0.5f, _k + 0.5f),
                                       _crns[6], mat);
        SpaceConverter::fromVoxToModel(point(_i - 0.5f, _j + 0.5f, _k + 0.5f),
                                       _crns[7], mat);
    }
    inline void get_corners_coords(int _i, int _j, int _k,
                                   const shared_ptr<Volume3DScalar>& _vol,
                                   point* _crns)
    {
        // const auto& mat = _vol->getVoxToModelMat();
        auto mat = trimesh::xform::identity();
        SpaceConverter::fromVoxToModel(point(_i + 0.5f, _j - 0.5f, _k - 0.5f),
                                       _crns[0], mat);
        SpaceConverter::fromVoxToModel(point(_i + 0.5f, _j + 0.5f, _k - 0.5f),
                                       _crns[1], mat);
        SpaceConverter::fromVoxToModel(point(_i - 0.5f, _j - 0.5f, _k - 0.5f),
                                       _crns[2], mat);
        SpaceConverter::fromVoxToModel(point(_i - 0.5f, _j + 0.5f, _k - 0.5f),
                                       _crns[3], mat);
        SpaceConverter::fromVoxToModel(point(_i + 0.5f, _j - 0.5f, _k + 0.5f),
                                       _crns[4], mat);
        SpaceConverter::fromVoxToModel(point(_i + 0.5f, _j + 0.5f, _k + 0.5f),
                                       _crns[5], mat);
        SpaceConverter::fromVoxToModel(point(_i - 0.5f, _j - 0.5f, _k + 0.5f),
                                       _crns[6], mat);
        SpaceConverter::fromVoxToModel(point(_i - 0.5f, _j + 0.5f, _k + 0.5f),
                                       _crns[7], mat);
    }
    // return the unique ids for corners at voxel ijk.
    inline void get_corners_ids(int _i, int _j, int _k,
                                const shared_ptr<Volume3DScalar>& _vol,
                                ivec3* _corner_ids)
    {
        _i++;
        _j++;
        _k++;
        _corner_ids[0] = ivec3(2 * _i + 1, 2 * _j - 1, 2 * _k - 1);
        _corner_ids[1] = ivec3(2 * _i + 1, 2 * _j + 1, 2 * _k - 1);
        _corner_ids[2] = ivec3(2 * _i - 1, 2 * _j - 1, 2 * _k - 1);
        _corner_ids[3] = ivec3(2 * _i - 1, 2 * _j + 1, 2 * _k - 1);
        _corner_ids[4] = ivec3(2 * _i + 1, 2 * _j - 1, 2 * _k + 1);
        _corner_ids[5] = ivec3(2 * _i + 1, 2 * _j + 1, 2 * _k + 1);
        _corner_ids[6] = ivec3(2 * _i - 1, 2 * _j - 1, 2 * _k + 1);
        _corner_ids[7] = ivec3(2 * _i - 1, 2 * _j + 1, 2 * _k + 1);
    }
    // return the "coordinates" of the 6 faces of voxel ijk.
    // each face is represented as 4-tuple of indices into the corner's indices
    // of the same voxel
    inline void get_faces_coords(int _i, int _j, int _k,
                                 const shared_ptr<Volume3DScalar>& _vol,
                                 ivec4* _faces)
    {
        _faces[0] = ivec4(0, 1, 3, 2);
        _faces[1] = ivec4(4, 5, 7, 6);
        _faces[2] = ivec4(0, 4, 5, 1);
        _faces[3] = ivec4(2, 6, 7, 3);
        _faces[4] = ivec4(4, 0, 2, 6);
        _faces[5] = ivec4(5, 1, 3, 7);
    }
    // return the unique ids for the 6 faces at voxel ijk
    inline void get_faces_ids(int _i, int _j, int _k,
                              const shared_ptr<Volume3DScalar>& _vol,
                              ivec3* _face_ids)
    {
        _i++;
        _j++;
        _k++;
        _face_ids[0] = ivec3(2 * _i, 2 * _j, 2 * _k - 1);
        _face_ids[1] = ivec3(2 * _i, 2 * _j, 2 * _k + 1);
        _face_ids[2] = ivec3(2 * _i + 1, 2 * _j, 2 * _k);
        _face_ids[3] = ivec3(2 * _i - 1, 2 * _j, 2 * _k);
        _face_ids[4] = ivec3(2 * _i, 2 * _j - 1, 2 * _k);
        _face_ids[5] = ivec3(2 * _i, 2 * _j + 1, 2 * _k);
    }
    // main init logic
    void init();
    // init neighbor voxel offset vectors (pointing from cur voxel to its
    // neighbor voxel)
    inline void set_offsets()
    {
        m_faceId_nbvoxOffset_map.resize(6);
        m_faceId_nbvoxOffset_map[0] = ivec3(-1, 0, 0);
        m_faceId_nbvoxOffset_map[1] = ivec3(1, 0, 0);
        m_faceId_nbvoxOffset_map[2] = ivec3(0, -1, 0);
        m_faceId_nbvoxOffset_map[3] = ivec3(0, 1, 0);
        m_faceId_nbvoxOffset_map[4] = ivec3(0, 0, -1);
        m_faceId_nbvoxOffset_map[5] = ivec3(0, 0, 1);
    }
    // init corners & face shared between a voxel and each of the 6 neighboring
    // voxel
    inline void set_elems_wrt_nbVox()
    {
        m_cornersWRTNbOffset.resize(6);
        m_cornersWRTNbOffset[0] = ivec4(2, 3, 7, 6);
        m_cornersWRTNbOffset[1] = ivec4(0, 1, 5, 4);
        m_cornersWRTNbOffset[2] = ivec4(4, 6, 2, 0);
        m_cornersWRTNbOffset[3] = ivec4(1, 3, 7, 5);
        m_cornersWRTNbOffset[4] = ivec4(0, 2, 3, 1);
        m_cornersWRTNbOffset[5] = ivec4(5, 7, 6, 4);
        m_faceWRTNbOffset.resize(6);
        m_faceWRTNbOffset[0] = 3;
        m_faceWRTNbOffset[1] = 2;
        m_faceWRTNbOffset[2] = 4;
        m_faceWRTNbOffset[3] = 5;
        m_faceWRTNbOffset[4] = 0;
        m_faceWRTNbOffset[5] = 1;
    }
    // init num of components per occupancy configuration of a cube's 8 corners
    // @note: should only be called by major init() entry
    void set_cube_cmpnts_num();
    // init for each edge of a voxel cube the offsets of the 4 voxels that share
    // the edge w.r.t. the cube center
    // @note: should only be called by major init() entry
    void set_nbvox_offsets_for_edges();
    // is the corner inside the volume
    inline bool is_corner_inside(const ivec3& _vox, int _c_id,
                                 const shared_ptr<Volume3DScalar>& _vol) const
    {
        point vox(_vox[0], _vox[1], _vox[2]);
        return SpaceConverter::allNbVoxelsInside(
            vox + SpaceConverter::vox_corner_offset[_c_id], _vol);
    }
    // is the edge of the voxel inside the volume
    inline bool is_edge_inside(const ivec3& _vox, int _e_id,
                               const shared_ptr<Volume3DScalar>& _vol) const
    {
        const auto& nb_vox_offset = m_nbVoxOffsetsForEdge[_e_id];
        int inside_cnt = 0;
        for (const auto& o : nb_vox_offset)
        {
            inside_cnt += SpaceConverter::get_occupancy_at_vox(_vox + o, _vol);
        }
        return inside_cnt == nb_vox_offset.size();
    }
    // is the face of the voxel inside the volume
    inline bool is_face_inside(const ivec3& _vox, int _f_id,
                               const shared_ptr<Volume3DScalar>& _vol) const
    {
        return SpaceConverter::get_occupancy_at_vox(
                   _vox + m_faceId_nbvoxOffset_map[_f_id], _vol) == 1;
    }
    // return whether two voxels have different values
    // out-of-bound voxel is assumed to take on value 0
    inline bool differ_in_values(const ivec3& _v1, const ivec3& _v2,
                                 const shared_ptr<Volume3DScalar>& _vol)
    {
        return SpaceConverter::get_occupancy_at_vox(_v1, _vol) !=
               SpaceConverter::get_occupancy_at_vox(_v2, _vol);
    }
    // triangulate the quad f and put resulting faces to _faces
    void triangulate_quad(const ivec4& _f, vector<uTriFace>& _faces);
    // save the voxel coords to the thin-layer set
    void save_to_thin_inside_layer(const ivec3& _vox,
                                   const shared_ptr<Volume3DScalar>& _vol);
    // find out how many connected components are there in the volume
    int compute_conn_cmpnt_num(const shared_ptr<Volume3DScalar>& _vol);

    /*data members*/
    // store the inside voxels' coords that are closest to the boundary
    // it's a thin layer of voxels, just one-voxel-thick
    set<ivec3> m_inside_voxs_thin;
    // num of conn. components
    int m_n_conn_cmpnts;
    // stores all vertices here
    vector<ivec3> m_vts;
    // stores all triangle faces here
    vector<uTriFace> m_tris;
    // pointer to the input volume rep
    std::shared_ptr<Volume3DScalar> m_vol;
    /*constants*/
    static bool consts_initialized;
    // how many surface components are there
    // given a standard marching-cube-8-corner-occupancy-bit-mask
    static int m_nCmpnts_per_cube[256];
    // mapping: face id -> the offset to get to the other voxel sharing the face
    static vector<ivec3> m_faceId_nbvoxOffset_map;
    //  the 4 corners' indices of that face
    static vector<ivec4> m_cornersWRTNbOffset;
    // mapping: index {0,..,5} -> the face id
    static vector<int> m_faceWRTNbOffset;
    // mapping: edge id -> 4 nb voxels that share the edge
    static vector<array<ivec3, 4>> m_nbVoxOffsetsForEdge;

protected: /*utility classes*/
    /* A walker that collects the centers of all LeafNodes of an octree */
    struct OctreeLeafCenterCollector : public OctreeVolume::Walker
    {
    protected:
        // the list of leaves' center
        vector<point>* m_centers;

    public:
        OctreeLeafCenterCollector() {}
        OctreeLeafCenterCollector(vector<point>& _centers)
            : m_centers(&_centers)
        {
        }
        ~OctreeLeafCenterCollector() {}
        // initialize the walker to the proper starting state
        void init();
        // override walker's default behavior
        void leaf(OctreeVolume* _tree, OctreeVolume::OctreeNode* _node,
                  const ivec3& _off, int _len);
    };
    /* A walker that collects boundary elements of all LeafNodes of an octree */
    struct OctreeBoundaryCollector : public OctreeVolume::Walker
    {
    protected:
    public:
    };
};

#endif