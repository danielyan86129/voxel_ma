#include <unordered_map>

#include <cassert>
#include <iostream>
#include <queue>
#include <set>

#include <isosurface/Cube.h> // Tao's cycle table

#include <voxelcore/densevolume.h>
#include <voxelcore/geomalgo.h>
#include <voxelcore/octree.h>
#include <voxelcore/reterrcode.h>
#include <voxelcore/surfacing.h>

using std::set;
using std::unordered_map;

/*
** static data members initialized here
*/
bool Surfacer::consts_initialized = false;
int Surfacer::m_nCmpnts_per_cube[] = {};
vector<ivec3> Surfacer::m_faceId_nbvoxOffset_map;
vector<ivec4> Surfacer::m_cornersWRTNbOffset;
vector<int> Surfacer::m_faceWRTNbOffset;
vector<array<ivec3, 4>> Surfacer::m_nbVoxOffsetsForEdge;

Surfacer::Surfacer() { init(); }

Surfacer::Surfacer(std::shared_ptr<Volume3DScalar> _vol)
{
    init();
    m_vol = _vol;
}

Surfacer::~Surfacer() {}

struct ivec3_less
{
    bool operator()(const ivec3& _a, const ivec3& _b) const
    {
        return _a[0] != _b[0]
                   ? _a[0] < _b[0]
                   : _a[1] != _b[1] ? _a[1] < _b[1]
                                    : _a[2] != _b[2] ? _a[2] < _b[2] : false;
    }
};
SurfacerErrCode
Surfacer::extractBoundaryVtsAndFaces(const shared_ptr<Volume3DScalar>& _vol,
                                     vector<point>& _vts, vector<ivec4>& _faces,
                                     bool _need_euler)
{
    if (!_vol)
        return SurfacerErrCode::EMPTY_VOL;
    m_inside_voxs_thin.clear();
    _vts.clear();
    _faces.clear();

    unordered_map<ivec3, int, ivec3Hash> access_elem_index;

    ivec3 cur_vox;
    ivec3 nb_vox;
    // ivec3 cnrs_coords[ 8 ];
    point cnrs_coords[8];
    ivec3 cnrs_ids[8];
    ivec4 faces_coords[6];
    ivec3 faces_ids[6];

    // extract the boundary bits around the given voxel
    auto process_voxel = [&](int i, int j, int k) {
        cur_vox = ivec3(i, j, k);
        get_corners_coords(i, j, k, _vol, cnrs_coords);
        get_corners_ids(i, j, k, _vol, cnrs_ids);
        get_faces_coords(i, j, k, _vol, faces_coords);
        get_faces_ids(i, j, k, _vol, faces_ids);
        for (auto o = 0; o < m_faceId_nbvoxOffset_map.size(); ++o)
        {
            nb_vox = cur_vox + m_faceId_nbvoxOffset_map[o];

            // // debug begin
            // auto cur_vox_val =
            //     SpaceConverter::get_occupancy_at_vox(cur_vox, _vol);
            // [[maybe_unused]] auto nb_vox_val =
            //     SpaceConverter::get_occupancy_at_vox(nb_vox, _vol);
            // if ( cur_vox_val == 0 /*&& nb_vox_val == 1 &&
            // 	(
            // 		nb_vox[ 0 ] == _vol->getSizeX() - 1 ||
            // 		nb_vox[ 1 ] == _vol->getSizeY() - 1 ||
            // 		nb_vox[ 2 ] == _vol->getSizeZ() - 1
            // 		)*/ )
            //     int stop = 1;
            // // debug end

            if (differ_in_values(cur_vox, nb_vox, _vol))
            {
                // found a pair of voxels sandwiching a boundary face
                // save the inside voxel coords
                if (_need_euler)
                {
                    SpaceConverter::get_occupancy_at_vox(cur_vox, _vol) == 1
                        ? save_to_thin_inside_layer(cur_vox, _vol)
                        : save_to_thin_inside_layer(nb_vox, _vol);
                }

                const auto& one_side_crnrs = m_cornersWRTNbOffset[o];
                const auto fi = m_faceWRTNbOffset[o];
                // for each corner between curvox and nbvox, test whether it
                // exists or not
                for (auto ii = 0; ii < one_side_crnrs.size(); ++ii)
                {
                    auto ci = one_side_crnrs[ii];
                    auto find_it = access_elem_index.find(cnrs_ids[ci]);
                    if (find_it == access_elem_index.end())
                    {
                        // this corner wasn't created before. add it to
                        // vts_list.
                        const auto& c = cnrs_coords[ci];
                        _vts.push_back(c);
                        //_vts.push_back( ivec3( c[ 0 ], c[ 1 ], c[ 2 ] ) );
                        access_elem_index[cnrs_ids[ci]] =
                            (int)(_vts.size() - 1);
                    }
                }
                // test whether the striding face exists or not
                if (access_elem_index.find(faces_ids[fi]) ==
                    access_elem_index.end())
                {
                    // this face wasn't there before. create it.
                    auto f = faces_coords[fi];
                    auto cid = cnrs_ids[f[0]];
                    auto find_it = access_elem_index.find(cid);
                    assert(find_it != access_elem_index.end());
                    f[0] = find_it->second;
                    cid = cnrs_ids[f[1]];
                    find_it = access_elem_index.find(cid);
                    assert(find_it != access_elem_index.end());
                    f[1] = find_it->second;
                    cid = cnrs_ids[f[2]];
                    find_it = access_elem_index.find(cid);
                    assert(find_it != access_elem_index.end());
                    f[2] = find_it->second;
                    cid = cnrs_ids[f[3]];
                    find_it = access_elem_index.find(cid);
                    assert(find_it != access_elem_index.end());
                    f[3] = find_it->second;
                    _faces.push_back(f);
                    access_elem_index[faces_ids[fi]] = (int)(_faces.size() - 1);
                }
            }
        }
    };

    if (std::dynamic_pointer_cast<DenseVolume>(_vol))
    { // process every voxel since we are dealing with dense volume
        for (auto i = 0; i < _vol->getSizeX(); ++i)
        {
            for (auto j = 0; j < _vol->getSizeY(); ++j)
            {
                for (auto k = 0; k < _vol->getSizeZ(); ++k)
                {
                    process_voxel(i, j, k);
                }
            }
        }
    }
    else if (auto vol_ptr = std::dynamic_pointer_cast<OctreeVolume>(_vol))
    { // process only boundary voxels since we are dealing
      // with octree
        vector<ivec3> boundary_voxels;
        cout << "getting boundary voxels..." << endl;
        vol_ptr->getBoundaryVoxels(boundary_voxels);
        cout << "# boundary voxels obtained: " << boundary_voxels.size()
             << endl;
        for (const auto& v : boundary_voxels)
        {
            process_voxel(v[0], v[1], v[2]);
        }
        cout << "all boundary bits extracted." << endl;
        boundary_voxels.clear();
        boundary_voxels.shrink_to_fit();
    }
    else
    { // volume type not supported
        return SurfacerErrCode::INVALID_VOL_TYPE;
    }

    access_elem_index.clear();

    // compute the num of connected components
    // using the just collected thin-inside-voxels
    if (_need_euler)
        m_n_conn_cmpnts = compute_conn_cmpnt_num(_vol);

    return SurfacerErrCode::SUCCESS;
}

SurfacerErrCode Surfacer::extractBoundaryVtsAndFaces(
    const shared_ptr<Volume3DScalar>& _vol, vector<point>& _vts,
    vector<uTriFace>& _faces, eulerchar* _euler)
{
    vector<ivec4> quads;
    auto errcode =
        extractBoundaryVtsAndFaces(_vol, _vts, quads, _euler != nullptr);
    if (errcode != SurfacerErrCode::SUCCESS)
        return errcode;

    if (_euler)
    {
        computeEulerChar(_vts, quads, _vol, *_euler);
    }

    _faces.clear();
    for (auto& f : quads)
    {
        triangulate_quad(f, _faces);
    }

    return SurfacerErrCode::SUCCESS;
}

SurfacerErrCode
Surfacer::extractBoundaryVts(const shared_ptr<Volume3DScalar>& _vol,
                             vector<point>& _vts)
{
    if (!_vol)
        return SurfacerErrCode::EMPTY_VOL;
    _vts.clear();

    unordered_map<ivec3, int, ivec3Hash> access_elem_index;

    ivec3 cur_vox;
    ivec3 nb_vox;
    // ivec3 cnrs_coords[ 8 ];
    point cnrs_coords[8];
    ivec3 cnrs_ids[8];

    // extract boundary vts around the given voxel
    auto process_voxel = [&](int _i, int _j, int _k) {
        cur_vox = ivec3(_i, _j, _k);
        get_corners_coords(_i, _j, _k, _vol, cnrs_coords);
        get_corners_ids(_i, _j, _k, _vol, cnrs_ids);
        for (auto o = 0; o < m_faceId_nbvoxOffset_map.size(); ++o)
        {
            nb_vox = cur_vox + m_faceId_nbvoxOffset_map[o];
            if (differ_in_values(cur_vox, nb_vox, _vol))
            {
                const auto& one_side_crnrs = m_cornersWRTNbOffset[o];
                [[maybe_unused]] const auto fi = m_faceWRTNbOffset[o];
                // for each corner between curvox and nbvox, test whether it
                // exists or not
                for (auto ii = 0; ii < one_side_crnrs.size(); ++ii)
                {
                    auto ci = one_side_crnrs[ii];
                    auto find_it = access_elem_index.find(cnrs_ids[ci]);
                    if (find_it == access_elem_index.end())
                    {
                        // this corner wasn't created before. add it to
                        // vts_list.
                        const auto& c = cnrs_coords[ci];
                        _vts.push_back(c);
                        //_vts.push_back( ivec3( c[ 0 ], c[ 1 ], c[ 2 ] ) );
                        access_elem_index[cnrs_ids[ci]] =
                            (int)(_vts.size() - 1);
                    }
                }
            }
        }
    };

    if (std::dynamic_pointer_cast<DenseVolume>(_vol))
    { // process every voxel since we are dealing with dense volume
        for (auto i = 0; i < _vol->getSizeX(); ++i)
        {
            for (auto j = 0; j < _vol->getSizeY(); ++j)
            {
                for (auto k = 0; k < _vol->getSizeZ(); ++k)
                {
                    process_voxel(i, j, k);
                }
            }
        }
    }
    else if (auto vol_ptr = std::dynamic_pointer_cast<OctreeVolume>(_vol))
    { // process only boundary voxels since we are dealing
      // with octree
        /*vector<ivec3> boundary_voxels;
        vol_ptr->getBoundaryVoxels(boundary_voxels);
        for ( const auto& v : boundary_voxels )
        {
                process_voxel( v[ 0 ], v[ 1 ], v[ 2 ] );
        }
        boundary_voxels.clear();
        boundary_voxels.shrink_to_fit();*/

        // fill the vertex list with a walker
        _vts.clear();
        OctreeLeafCenterCollector walker(_vts);
        vol_ptr->walkTree(&walker);
        cout << "# leaf centers collected: " << _vts.size() << endl;

        /// UNCOMMENT following to convert centers from voxel to model space
        ///  for now don't convert and leave the decision to the caller
        // point p;
        // for ( auto& v : _vts )
        //{
        //	//SpaceConverter::fromVoxToModel( vox_cnr, p,
        //_tree->getVoxToModelMat()
        //); 	SpaceConverter::fromVoxToModel( v, p, _vol->getVoxToModelMat()
        //); v = p;
        // }
    }
    else
    { // volume type not supported
        return SurfacerErrCode::INVALID_VOL_TYPE;
    }

    return SurfacerErrCode::SUCCESS;
}

int Surfacer::computeEulerChar(const vector<point>& _vts,
                               const vector<ivec4>& _faces,
                               const shared_ptr<Volume3DScalar>& _vol,
                               eulerchar& _euler) const
{
    int euler = 0;

    /* compute the euler char for the boundary part first */
    unordered_map<ivec2, int, ivec2Hash> edge_refCntByFaces;
    // predicate for deciding whether edge/vert is singular
    auto multiplicity_e = [&edge_refCntByFaces](const ivec2 _e) -> int {
        return edge_refCntByFaces.find(_e)->second == 2 ? 1 : 2;
    };
    auto multiplicity_v = [&](int _vi) -> int {
        /*if ( vert_refCntByFaces[ _vi ] != 6 ||
                vert_refCntByFaces[ _vi ] != 9 )
                return 1;*/

        // grab the 8 nb voxels values from _vol
        const auto& v = _vts[_vi];
        point p((float)v[0], (float)v[1], (float)v[2]);
        point q;
        ivec3 nb_voxels[8];
        SpaceConverter::fromModelToVox(
            p, q, /*_vol->getModeltoVoxMat()*/ trimesh::xform::identity());
        SpaceConverter::getNbVoxelCoords(q, nb_voxels);
        unsigned char mask = 0;
        for (int i = 0; i < 8; ++i)
        {
            mask |= (SpaceConverter::get_occupancy_at_vox(nb_voxels[i], _vol)
                     << (7 - i));
        }

        auto n_conn_cmpnt = m_nCmpnts_per_cube[mask];
        return n_conn_cmpnt;
    };

    for (const auto& f : _faces)
    {
        edge_refCntByFaces[util::makeEdge(f[0], f[1])]++;
        edge_refCntByFaces[util::makeEdge(f[1], f[2])]++;
        edge_refCntByFaces[util::makeEdge(f[2], f[3])]++;
        edge_refCntByFaces[util::makeEdge(f[3], f[0])]++;
    }

    // V E F
    int V = 0, E = 0;
    auto F = _faces.size();
    // compute singular vert and edge number
    int V_s = 0, E_s = 0;
    for (int vi = 0; vi < _vts.size(); ++vi)
    {
        int mult_v = multiplicity_v(vi);
        if (mult_v != 1)
            V_s++;
        V += mult_v;
    }
    for (const auto& e_ref_pr : edge_refCntByFaces)
    {
        int mult_e = multiplicity_e(e_ref_pr.first);
        if (mult_e > 1)
            E_s++;
        E += mult_e;
    }
    std::cout << "# singular V/E = " << V_s << "/" << E_s << std::endl;
    edge_refCntByFaces.clear();
    /* boundary part euler char computed */

    /*now compute the euler char of the interior part*/
    int V_int = 0, E_int = 0, F_int = 0;
    int T = 0; // the volume element (# cubes)
    for (int i = 0; i < _vol->getSizeX(); ++i)
    {
        for (int j = 0; j < _vol->getSizeY(); ++j)
        {
            for (int k = 0; k < _vol->getSizeZ(); ++k)
            {
                const auto& v = ivec3(i, j, k);
                if (SpaceConverter::get_occupancy_at_vox(v, _vol) == 1)
                {
                    // count this interior cube
                    T++;
                    // count interior corners
                    for (int ii = 0; ii < 8; ++ii)
                    {
                        V_int += (int)is_corner_inside(v, ii, _vol);
                    }
                    // count interior edges
                    for (int ii = 0; ii < 12; ++ii)
                    {
                        E_int += (int)is_edge_inside(v, ii, _vol);
                    }
                    // count interior faces
                    for (int ii = 0; ii < 6; ++ii)
                    {
                        F_int += (int)is_face_inside(v, ii, _vol);
                    }
                }
            }
        }
    } // done iterating thru all voxels
    // corners, edges, and faces are duplicated during counting.
    V_int /= 8;
    E_int /= 4;
    F_int /= 2;
    std::cout << "# interior V / E / F = " << V_int << " / " << E_int << " / "
              << F_int << std::endl;

    // now we need to figure out the last part of Euler Char: the num of
    // connected components
    auto C = m_n_conn_cmpnts;

    euler = (int)(V + V_int - E - E_int + F + F_int - T);

    _euler = eulerchar(V + V_int, E + E_int, F + F_int, T, C, euler);

    return euler;
}

/*
** helpers
*/

void Surfacer::init()
{
    if (!consts_initialized)
    {
        set_offsets();
        set_elems_wrt_nbVox();
        set_cube_cmpnts_num();
        set_nbvox_offsets_for_edges();
    }
    consts_initialized = true;
}

void Surfacer::set_cube_cmpnts_num()
{
    using std::queue;

    Cube cycl_tbl("cycle8.txt");
    vector<vector<int>> edge_edge_adj_tbl;
    queue<int> cmpnt_q;
    vector<bool> visited;
    for (unsigned char mask = 0; mask <= 255; ++mask)
    {
        // Quick test
        if (mask == 0 || mask == 255)
        {
            m_nCmpnts_per_cube[mask] = 0;
            if (mask == 255)
                break;
            continue;
        }

        // initially, each mask config corresponds to at least 1 surface
        // component and no edge is ref.ed by any cycle yet
        m_nCmpnts_per_cube[mask] = 0;
        edge_edge_adj_tbl.assign(12, vector<int>());

        auto cycl_list = cycl_tbl.getCycle(mask);
        int numcycls = cycl_list->getLength();
        Cycle* cycl;
        int cycl_len;
        int cur_e_id, pre_e_id, next_e_id;

        // first, loop through all cycles to count
        // for each edge how many times they have been ref.ed
        for (int i = 0; i < numcycls; ++i)
        {
            cycl = cycl_list->getFirst();
            cycl_len = cycl->getLength();
            for (int j = 0; j < cycl_len; ++j)
            {
                cur_e_id = cycl->getEdge(j);
                if (j > 0)
                {
                    // put prev edge to adjacency
                    pre_e_id = cycl->getEdge(j - 1);
                    edge_edge_adj_tbl[cur_e_id].push_back(pre_e_id);
                }
                if (j < cycl_len - 1)
                {
                    // put next edge to adjacency
                    next_e_id = cycl->getEdge(j + 1);
                    edge_edge_adj_tbl[cur_e_id].push_back(next_e_id);
                }
            }
            cycl_list->rotateLeft();
        }

        // second, trace out the connected components
        visited.assign(12, false);
        int ncc = 0;
        for (int i = 0; i < 12; ++i)
        {
            if (edge_edge_adj_tbl[i].empty() || visited[i])
                continue;

            // new component starts here
            ncc++;
            cmpnt_q.push(i);
            while (!cmpnt_q.empty())
            {
                int cur = cmpnt_q.front();
                cmpnt_q.pop();
                if (visited[cur])
                    continue;
                visited[cur] = true;
                // add cur's neighbors (from adj table) to queue
                const auto& nbs = edge_edge_adj_tbl[cur];
                for (auto ni : nbs)
                {
                    if (!visited[ni])
                    {
                        cmpnt_q.push(ni);
                    }
                }
            }
        }
        // now we have the # connected surface components for this mask config
        m_nCmpnts_per_cube[mask] = ncc;
    }
    // if the cube has any of the following config,
    // then the center of the cube is a "hole"
    m_nCmpnts_per_cube[(unsigned char)126] = 0;
    m_nCmpnts_per_cube[(unsigned char)189] = 0;
    m_nCmpnts_per_cube[(unsigned char)219] = 0;
    m_nCmpnts_per_cube[(unsigned char)231] = 0;
    std::cout << "cycle table set." << std::endl;
}

void Surfacer::set_nbvox_offsets_for_edges()
{
    m_nbVoxOffsetsForEdge.resize(12);
    // possible config of ijk. Each seq. specifies which one of ijk should stay
    // 0, followed by the other 2 as changing iterators
    ivec3 ijk_seq[3] = {ivec3(0, 1, 2), ivec3(1, 0, 2), ivec3(2, 1, 0)};
    ivec3 ijk;
    // two possible ranges for the other 2 to iterate over
    vector<ivec2> ranges(2);
    ranges[0] = ivec2(-1, 0);
    ranges[1] = ivec2(0, 1);
    // keep track of current edge id and which nb vox of it
    int ei = 0, ni = 0;
    for (const auto& s : ijk_seq)
    {
        int const_idx = s[0];
        int j = s[1], k = s[2];
        ijk[const_idx] = 0;
        for (const auto& r1 : ranges)
        {
            for (const auto& r2 : ranges)
            {
                ni = 0;
                for (ijk[j] = r1[0]; ijk[j] <= r1[1]; ++ijk[j])
                {
                    for (ijk[k] = r2[0]; ijk[k] <= r2[1]; ++ijk[k])
                    {
                        m_nbVoxOffsetsForEdge[ei][ni] = ijk;
                        ni++;
                    }
                }
                ei++;
            }
        }
    }
}

void Surfacer::triangulate_quad(const ivec4& _f, vector<uTriFace>& _faces)
{
    _faces.push_back(uTriFace(_f[0], _f[1], _f[2]));
    _faces.push_back(uTriFace(_f[2], _f[3], _f[0]));
}

void Surfacer::save_to_thin_inside_layer(const ivec3& _vox,
                                         const shared_ptr<Volume3DScalar>& _vol)
{
    // save all possible voxels among 27 local neighbors that are inside
    ivec3 nb;
    for (int i = -1; i <= 1; ++i)
    {
        for (int j = -1; j <= 1; ++j)
        {
            for (int k = -1; k <= 1; ++k)
            {
                nb = _vox + ivec3(i, j, k);
                if (SpaceConverter::get_occupancy_at_vox(_vox + ivec3(i, j, k),
                                                         _vol) == 1)
                {
                    m_inside_voxs_thin.insert(nb);
                }
            }
        }
    }
}

int Surfacer::compute_conn_cmpnt_num(const shared_ptr<Volume3DScalar>& _vol)
{
    using std::queue;

    int ncc = 0;
    set<ivec3> visited;
    queue<ivec3> cmpnt_q;

    // use 6-connected way to find connected components
    ivec3 nb;
    for (const auto& v : m_inside_voxs_thin)
    {
        // find an unvisited vox to start with
        if (visited.count(v) > 0)
            continue;
        // now start a new component
        ncc++;
        cmpnt_q.push(v);
        while (!cmpnt_q.empty())
        {
            auto cur_v = cmpnt_q.front();
            cmpnt_q.pop();
            if (visited.count(cur_v) > 0)
                continue;
            // else, set it as visited, and expand the unvisited & inside
            // neighbors of cur_v into the queue
            visited.insert(cur_v);
            for (const auto& o : m_faceId_nbvoxOffset_map)
            {
                nb = cur_v + o;
                // if this nb is a inside vox & not yet visited
                // then visit it
                if (m_inside_voxs_thin.count(nb) > 0 && visited.count(nb) == 0)
                {
                    cmpnt_q.push(nb);
                }
            }
        }
    }

    return ncc;
}

/*
private helper classes
*/

void Surfacer::OctreeLeafCenterCollector::init() { m_centers->clear(); }

void Surfacer::OctreeLeafCenterCollector::leaf(OctreeVolume* _tree,
                                               OctreeVolume::OctreeNode* _node,
                                               const ivec3& _off, int _len)
{
    /*assert( _node->getType() == OctreeVolume::NodeType::LEAF &&
            "Strange: leaf() called on non-leaf node..." );*/
    if (_node->getType() != OctreeVolume::LEAF)
        return;

    auto leafvalue = dynamic_cast<OctreeVolume::LeafNode*>(_node)->getValues();
    if (leafvalue == (unsigned char)0 || leafvalue == (unsigned char)0xFF)
        return;
    // the center of a node is really just
    // the (1, 1, 1) corner of the lowest voxel of the node
    auto vox_cnr = SpaceConverter::getVoxelCorner(_off, 1, 1, 1);
    ivec3 p;
    // SpaceConverter::fromVoxToModel( vox_cnr, p, _tree->getVoxToModelMat() );
    m_centers->push_back(vox_cnr);
}
