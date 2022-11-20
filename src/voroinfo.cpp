#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>

#include <ANN/ANN.h>
#include <trimesh/TriMesh_algo.h>

#include <voxelcore/geomalgo.h>
#include <voxelcore/voroinfo.h>

namespace voxelvoro
{

using std::cout;
using std::endl;
using std::ifstream;
using std::numeric_limits;
using std::ofstream;
using std::queue;
using std::set;
using std::string;
using trimesh::vec;

/********************************************** */
/* class methods implementation
/********************************************** */
VoroInfo::VoroInfo()
{
    // m_r_valid = false;
    invalidate_geometric_states();
}

VoroInfo::~VoroInfo() {}

bool VoroInfo::loadFromTetgenFiles(const char* _file_basename,
                                   shared_ptr<Volume3DScalar> _vol)
{
    timer t_IO, t_infer_sites;
    // at least: .v.node, .v.edge, .v.face need to be available
    ifstream in_node(string(_file_basename) + ".v.node");
    ifstream in_edge(string(_file_basename) + ".v.edge");
    ifstream in_face(string(_file_basename) + ".v.face");
    // ifstream in_cells( string( _file_basename ) + ".v.cell" );
    if (!in_node || !in_edge || !in_face)
    {
        cout << "Error: some input files (.v.node/edge/face/cell) are missing."
             << endl;
        return false;
    }

    invalidate_geometric_states();

    unordered_map<int, int>* old_new_v_map =
        _vol ? new unordered_map<int, int>() : nullptr;
    unordered_map<int, int>* old_new_e_map =
        _vol ? new unordered_map<int, int>() : nullptr;
    t_IO.start();
    load_voro_vts(in_node, _vol, old_new_v_map);
    cout << ".v.node file processed." << endl;
    vector<ivec2> edges;
    load_voro_edges(in_edge, _vol, old_new_v_map, old_new_e_map);
    cout << ".v.edge file processed." << endl;
    load_voro_faces(in_face, _vol, old_new_e_map);
    cout << ".v.face file processed." << endl;
    if (old_new_v_map)
        delete old_new_v_map;
    if (old_new_e_map)
        delete old_new_e_map;
    t_IO.stop();

    if (_vol)
        set_only_inside(true);
    /*
    ** uncomment the following if correspondence between voro cell index and
    *site index is unknown
    */
    // t_infer_sites.start();
    // infer_sites_from_cell_file( in_cells );
    // t_infer_sites.stop();
    // cout << ".v.cell file processed." << endl;
    /*
    ** The correspondence is known at least for voro output from tetgen-1.5
    *(i.g. cell i is dual to site i)
    ** so we don't need cell information to compute site-related info
    */
    cout << "computing info related to sites..." << endl;
    computeInfoRelatedtoSites();
    cout << "time(I/O) -> load voro from tetgen files: "
         << t_IO.elapseMilli().count() << " ms" << endl;
    // cout << "time -> infer sites from cell: " <<
    // t_infer_sites.elapseMilli().count() << " ms" << endl;

    /* TODO: make sure this always works
    ** Here we optimize away this finalize(), as it is redundant:
    ** for now we are sure the operations entail will only use basic geometry
    *info
    ** mostly, this voro will be simplified by merging operation, so most info
    ** computed by finalize() is gone anyway
    */
    // m_geom.finalize();

    return true;
} // loadFromTetgenFiles()

bool VoroInfo::loadFromTetgenFiles(const tetgenio& _tetio,
                                   shared_ptr<Volume3DScalar> _vol)
{
    timer t;
    invalidate_geometric_states();

    t.start();

    unordered_map<int, int> old_new_v_map;
    unordered_map<int, int> old_new_e_map;
    // load voro verts
    m_geom.clearVertList();
    this->m_is_finite_v.clear();
    const REAL* vpts = _tetio.vpointlist;
    for (size_t i = 0; i < _tetio.numberofvpoints; ++i)
    {
        point v(vpts[i * 3], vpts[i * 3 + 1], vpts[i * 3 + 2]);
        if (_vol && this->tagVert(v, _vol) == false)
            continue;
        if (_vol)
            old_new_v_map[i] = m_geom.numVts();
        m_geom.appendVert(v);
        if (!_vol)
            m_is_finite_v.push_back(true);
        m_bbox += v;
    }
    cout << "Voro nodes loaded." << endl;

    // load voro edges
    m_geom.clearEdgeList();
    m_is_finite_e.clear();
    for (size_t eid = 0; eid < _tetio.numberofvedges; ++eid)
    {
        const auto& ve = _tetio.vedgelist[eid];
        ivec2 e(ve.v1, ve.v2);
        vec dir;        // dir for infinite edge
        if (e[1] == -1) // infinite edge
            dir = vec(ve.vnormal);

        if (_vol) // only keep inside part if vol is present
        {
            if (!old_new_v_map.count(e[0]) || !old_new_v_map.count(e[1]))
                continue;
            e[0] = old_new_v_map[e[0]];
            e[1] = old_new_v_map[e[1]];
            old_new_e_map[eid] = m_geom.numEdges();
        }
        if (!_vol)
            if (e[1] == -1)
            {
                m_is_finite_e.push_back(false);
                // create an infinite vertex for this infinite edge to refer to
                point inf_v =
                    geom().getVert(e[0]) - dir * m_bbox.radius() * 0.01f;
                m_geom.appendVert(inf_v);
                m_is_finite_v.push_back(false);
                e[1] = m_geom.numVts() - 1;
            }
            else
                m_is_finite_e.push_back(true);
        m_geom.appendEdge(e);
    }

    // load voro faces
    m_geom.clearFaceList();
    m_face_sites.clear();
    m_is_finite_f.clear();
    vector<int> f_erep;
    vector<int> f_vts;
    int s1, s2, nsides, eid;
    for (size_t i = 0; i < _tetio.numberofvfacets; ++i)
    {
        const auto& vf = _tetio.vfacetlist[i];
        nsides = vf.elist[0];
        // in tetgen output, cell-id == site-id
        s1 = vf.c1;
        s2 = vf.c2;
        // read in n-edges
        f_erep.clear();
        for (auto j = 1; j <= nsides; ++j)
        {
            f_erep.push_back(vf.elist[j]);
        }
        // only keep inside face if vol present
        if (_vol)
            for (auto j = 0; j < f_erep.size(); ++j)
            {
                eid = f_erep[j];
                if (!old_new_e_map.count(eid))
                    goto NEXT_FACE;
                else
                    eid = old_new_e_map[eid];
                f_erep[j] = eid;
            }
        else // else, need to deal with infinite side
        {
            if (f_erep.back() == -1)
            {
                // use adjacent 2 half-infinite edges to create this new edge
                auto e1 = geom().getEdge(f_erep.front());
                auto e2 = geom().getEdge(*(f_erep.end() - 2));
                assert(!isVertexFinite(e1[1]) && !isVertexFinite(e2[1]));
                auto new_e = util::makeEdge(e1[1], e2[1]);
                m_geom.appendEdge(new_e);
                m_is_finite_e.push_back(false);
                f_erep.back() = geom().numEdges() - 1;
                m_is_finite_f.push_back(false);
            }
            else
                m_is_finite_f.push_back(true);
        }
        f_vts.resize(nsides);
        trace_face(f_erep, m_geom.m_edges, f_vts);
        m_geom.appendFace(f_vts);
        // save sites for this face
        m_face_sites.push_back(ivec2(s1, s2));

    NEXT_FACE: // refer to this label to quickly skip this face
               ;
    }
    m_face_sites.shrink_to_fit();
    m_is_finite_f.shrink_to_fit();
    if (_vol)
        set_only_inside(true);
    cout << "computing info related to sites..." << endl;
    computeInfoRelatedtoSites();

    t.stop();

    cout << "time -> load voro from tetgen: " << t.elapseMilli().count()
         << " ms" << endl;
    return true;
}

void VoroInfo::setInfo(const vector<point>& _vts, const vector<ivec2>& _edges,
                       const vector<vector<int>> _faces,
                       const vector<bool>& _vert_valid_tag,
                       const vector<bool>& _is_finite_v,
                       const vector<bool>& _is_finite_e,
                       const vector<bool>& _is_finite_f)
{
    m_geom.clear();
    for (const auto& v : _vts)
    {
        m_geom.appendVert(v);
    }
    for (const auto& e : _edges)
    {
        m_geom.appendEdge(e);
    }
    for (const auto& f : _faces)
    {
        m_geom.appendFace(f);
    }
    m_geom.finalize();
    m_vts_valid = _vert_valid_tag;
    m_is_finite_v = _is_finite_v;
    m_is_finite_e = _is_finite_e;
    m_is_finite_f = _is_finite_f;
}

void VoroInfo::setSitesPositions(const vector<point>& _sites_pos)
{
    // save the sites' position
    m_site_positions = _sites_pos;
}

void VoroInfo::computeInfoRelatedtoSites()
{
    m_face_sites_valid = true;

    // 2: estimate radius function
    // from the sites indices of each face derive the sites for each vertex
    // for each face f, pick a random site s, then,
    //     for each v of f, compute dist(v, s), assign that to r-of-v
    m_r_per_v.resize(m_geom.numVts(), 0.0f);
    vector<int> vts_f;
    float eps = getInsidePartSize() * 1.0e-5f;
    int inconst_cnt = 0;
    for (int fi = 0; fi < m_face_sites.size(); ++fi)
    {
        m_geom.getFaceVRep(fi, vts_f);
        int site_i = m_face_sites[fi][0]; // just pick the first site
        auto site_p = m_site_positions[site_i];
        for (auto vi : vts_f)
        {
            auto v_p = m_geom.getVert(vi);
            auto r = trimesh::dist(site_p, v_p);
            auto old_r = m_r_per_v[vi];
            if (old_r > 0.0f && !util::is_equal(r, old_r, eps))
            {
                /*cout << "Potential bug!! inconsistent radius. (for voro vert:
                " << vi
                << ", "
                << "old_r: " << old_r << ", new_r: " << r << endl;*/
                inconst_cnt++;
            }
            m_r_per_v[vi] = r;
        }
    }
    if (inconst_cnt)
        cout << "Potential bug!! inconsistent radius (at " << inconst_cnt
             << " voro vts)" << endl;

    cout << "radii estimated." << endl;
    m_r_valid = true;
}
void VoroInfo::computeInfoRelatedtoSites(const vector<point>& _cell_cents)
{
    int n_pts = m_site_positions.size();
    float sq_range = m_bbox.radius() * m_bbox.radius();
    ANNpointArray data_pts;
    ANNpoint query_p;
    ANNdistArray dists;
    ANNidxArray nn_idx;
    ANNkd_tree* kdtree;
    query_p = annAllocPt(3);
    data_pts = annAllocPts(n_pts, 3);
    for (size_t i = 0; i < n_pts; ++i)
    {
        const auto& v = m_site_positions[i];
        data_pts[i][0] = v[0];
        data_pts[i][1] = v[1];
        data_pts[i][2] = v[2];
    }
    kdtree = new ANNkd_tree(data_pts, n_pts, 3);
    nn_idx = new ANNidx[1];
    dists = new ANNdist[1];
    cout << "done: " << n_pts << " site positions inserted to kdtree." << endl;

    // 1. find corresponding sites for each face using input list of cell
    // centroids (right now m_face_sites contain index into this list)
    vector<int> cell_site_id(_cell_cents.size(), -1);
    for (auto ci = 0; ci < _cell_cents.size(); ++ci)
    {
        const auto& cent = _cell_cents[ci];
        query_p[0] = cent[0];
        query_p[1] = cent[1];
        query_p[2] = cent[2];
        /*auto nn_within_r = kdtree->annkFRSearch( query_p, sq_range, 1, nn_idx,
        nullptr, 0.0 ); if ( nn_within_r == 0 ) cout << "WARNING: cell centroid
        "<<ci<<" has zero closest site. Impossible!" << endl; else
        cell_site_id[ ci ] = nn_idx[ 0 ];*/
        kdtree->annkSearch(query_p, 1, nn_idx, dists, 0.0);
        cell_site_id[ci] = nn_idx[0];
    }
    // now collect sites into each face
    for (auto& sts : m_face_sites)
    {
        sts[0] = cell_site_id[sts[0]];
        sts[1] = cell_site_id[sts[1]];
    }
    m_face_sites_valid = true;

    // 2: estimate radius function
    m_r_per_v.resize(m_geom.numVts(), 0.0f);

    //// Solution 1: find nearest boundary point for each ma point
    //// their distance is the radius for that ma point
    //// by default each vertex associated with radius of 0.0
    // for ( int vi = 0; vi < m_geom.numVts(); ++vi )
    //{
    //	if ( /*isVertexFinite( vi )*/true )
    //	{
    //		const auto& p = m_geom.getVert( vi );
    //		query_p[ 0 ] = p[ 0 ];
    //		query_p[ 1 ] = p[ 1 ];
    //		query_p[ 2 ] = p[ 2 ];
    //		/*
    //		auto nn_within_r = kdtree->annkFRSearch( query_p, sq_range, 1,
    // nn_idx, nullptr, 0.0 ); 		if ( nn_within_r == 0 )
    // cout
    // << "WARNING: vert " << vi << ", zero NN found!" << endl;
    // else m_r_per_v[ vi ] = trimesh::dist( geom().getVert( nn_idx[ 0 ] ), p
    // );*/ kdtree->annkSearch( query_p, 1, nn_idx, dists, 0.0 );
    // m_r_per_v[ vi ] = trimesh::dist( geom().getVert( nn_idx[ 0 ] ), p );
    //	}
    // }

    // Solution 2: from the sites indices of each face derive the sites for each
    // vertex for each face f, pick a random site s, then,
    //     for each v of f, compute dist(v, s), assign that to r-of-v
    vector<int> vts_f;
    float eps = getInsidePartSize() * 1.0e-5f;
    int inconst_cnt = 0;
    for (int fi = 0; fi < m_face_sites.size(); ++fi)
    {
        m_geom.getFaceVRep(fi, vts_f);
        int site_i = m_face_sites[fi][0]; // just pick the first site
        auto site_p = m_site_positions[site_i];
        for (auto vi : vts_f)
        {
            auto v_p = m_geom.getVert(vi);
            auto r = trimesh::dist(site_p, v_p);
            auto old_r = m_r_per_v[vi];
            if (old_r > 0.0f && !util::is_equal(r, old_r, eps))
            {
                /*cout << "Potential bug!! inconsistent radius. (for voro vert:
                   " << vi
                   << ", "
                        << "old_r: " << old_r << ", new_r: " << r << endl;*/
                inconst_cnt++;
            }
            m_r_per_v[vi] = r;
        }
    }
    if (inconst_cnt)
        cout << "Potential bug!! inconsistent radius (at " << inconst_cnt
             << " voro vts)" << endl;

    cout << "radii estimated." << endl;
    m_r_valid = true;

    delete[] nn_idx;
    delete[] dists;
    delete kdtree;
    annClose();
}

const vector<point>& VoroInfo::getSitesPosition() const
{
    return m_site_positions;
}

bool VoroInfo::radiiValid() const { return m_r_valid; }

const vector<float>& VoroInfo::getRadii() const { return m_r_per_v; }

bool VoroInfo::tagVert(const point& _p,
                       const shared_ptr<Volume3DScalar>& _vol) const
{
    ivec3 vox; // the voxel containing a given 3d point
    SpaceConverter::fromModelToVox(
        _p, vox, /*_vol->getModeltoVoxMat()*/ trimesh::xform::identity());
    return SpaceConverter::voxTaggedAsInside(vox, _vol);
}

bool VoroInfo::tagVtsUsingUniformVol(const shared_ptr<Volume3DScalar>& _vol)
{
    m_vts_valid.resize(m_geom.numVts(), false);
    ivec3 vox; // the voxel containing a given 3d point
    for (size_t i = 0; i < m_geom.numVts(); ++i)
    {
        // find the voxel that point lies in
        const auto& pt = m_geom.getVert(i);
        SpaceConverter::fromModelToVox(
            pt, vox, /*_vol->getModeltoVoxMat()*/ trimesh::xform::identity());

        // the point will be invalid if the voxel is an "outside" voxel
        this->m_vts_valid[i] = SpaceConverter::voxTaggedAsInside(vox, _vol);
        //_vol->getDataAt( vox[ 0 ], vox[ 1 ], vox[ 2 ] ) >= 0.0;
    }
    m_edge_valid.resize(m_geom.numEdges(), false);
    for (auto i = 0; i < m_geom.numEdges(); ++i)
    {
        this->m_edge_valid[i] = computeEdgeValidity(i);
    }
    m_face_valid.resize(m_geom.numFaces(), false);
    for (auto i = 0; i < m_geom.numFaces(); ++i)
    {
        m_face_valid[i] = computeFaceValidity(i);
    }

    // update bbox
    m_bbox.clear();
    for (size_t i = 0; i < m_geom.numVts(); ++i)
    {
        if (isVertexValid(i))
        {
            auto& pt = m_geom.getVert(i);
            m_bbox += pt;
        }
    }

    return true;
}

float VoroInfo::getInsidePartSize() const { return m_bbox.size().max(); }

void VoroInfo::getValidVts(vector<int>& _vts_indices) const
{
    _vts_indices.clear();
    for (size_t vi = 0; vi < m_geom.numVts(); ++vi)
    {
        if (isVertexValid(vi))
            _vts_indices.push_back(vi);
    }
}

void VoroInfo::getValidEdges(vector<int>& _edge_indices) const
{
    _edge_indices.clear();
    for (size_t ei = 0; ei < m_geom.numEdges(); ++ei)
        if (isEdgeValid(ei))
            _edge_indices.push_back(ei);
}

void VoroInfo::getValidFaces(vector<int>& _face_indices) const
{
    _face_indices.clear();
    // for each face f, if one of its vts is invalid, f is not valid
    vector<int> f_v_rep;
    for (size_t i = 0; i < m_geom.numFaces(); ++i)
        if (isFaceValid(i))
            _face_indices.push_back(i);
}

void VoroInfo::getFiniteFaces(vector<int>& _face_indices) const
{
    _face_indices.clear();
    _face_indices.reserve(m_geom.numFaces());
    for (int fi = 0; fi < m_geom.numFaces(); ++fi)
    {
        if (isFaceFinite(fi))
            _face_indices.push_back(fi);
    }
    _face_indices.shrink_to_fit();
}

void VoroInfo::getAllFaces(vector<int>& _face_indices) const
{
    _face_indices.clear();
    _face_indices.reserve(m_geom.numFaces());
    for (int i = 0; i < m_geom.numFaces(); ++i)
        _face_indices.push_back(i);
}

void VoroInfo::getVts(const vector<int>& _vts_indices,
                      vector<point>& _vts) const
{
    m_geom.getVertCoords(_vts_indices, _vts);
}

void VoroInfo::getEdges(const vector<int>& _edges_indices,
                        vector<ivec2>& _edges) const
{
    _edges.reserve(_edges_indices.size());
    for (auto it = _edges_indices.begin(); it != _edges_indices.end(); ++it)
        _edges.push_back(m_geom.getEdge(*it));
}

void VoroInfo::getFaces(const vector<int>& _faces_indices,
                        vector<uTriFace>& _tri_faces) const
{
    vector<int> f_vrep;
    vector<uTriFace> tri_of_f;
    for (auto it = _faces_indices.begin(); it != _faces_indices.end(); ++it)
    {
        m_geom.getFaceVRep(*it, f_vrep);
        util::simpleTriangulate(f_vrep, tri_of_f);
        _tri_faces.insert(_tri_faces.end(), tri_of_f.begin(), tri_of_f.end());
    }
}

void VoroInfo::getPartInside(vector<int>& _vts_indices, vector<ivec2>& _edges,
                             vector<uTriFace>& _tri_faces) const
{
    // prepare vertices:
    // - compact away invalid vertices
    // - maintain a mapping: old vertex index -> new index in the list of
    // remaining vts
    vector<size_t> vert_map(m_geom.numVts());
    for (size_t vi = 0; vi < m_geom.numVts(); ++vi)
    {
        if (isVertexValid(vi))
        {
            vert_map[vi] = _vts_indices.size();
            _vts_indices.push_back(vi);
        }
    }

    // prepare edges and faces (triangulated):
    // - compact away invalid edges & faces,
    // - let edges and faces refer to newly remained vertices
    set<ivec2> total_tri_edges;
    vector<int> vts_of_f;
    vector<uTriFace> tri_faces_of_poly;
    ivec2 tri_edges_of_poly[3];
    uTriFace tcpy;
    for (size_t fi = 0; fi < m_geom.numFaces(); ++fi)
    {
        m_geom.getFaceVRep(fi, vts_of_f);
        util::simpleTriangulate(vts_of_f, tri_faces_of_poly);
        for (const auto& t : tri_faces_of_poly)
        {
            if (computeFaceValidity(t) == true)
            {
                // update index of triangle vertices
                tcpy[0] = vert_map[t[0]];
                tcpy[1] = vert_map[t[1]];
                tcpy[2] = vert_map[t[2]];
                // add the updated tri face
                _tri_faces.push_back(tcpy);
                // add the tri-edges
                util::makeEdgesFromTri(tcpy, tri_edges_of_poly);
                for (const auto& e : tri_edges_of_poly)
                    total_tri_edges.insert(e);
            }
        }
    }
    _edges.insert(_edges.end(), total_tri_edges.begin(), total_tri_edges.end());
    total_tri_edges.clear();
}

// TODO: optimize away redundant copy of msure, v/e/f, etc.
void VoroInfo::extractInsideWithMeasure(
    MeasureForMA::meassuretype _mssure_tp, vector<point>& _output_vts,
    vector<ivec2>& _output_edges, vector<uTriFace>& _output_tris,
    vector<int>& _from_fi, vector<float>& _vts_msure,
    vector<float>& _edges_msure, vector<float>& _faces_msure) const
{
    timer t_msure, t_get_inside;

    // grab indices to inside elements, i.e. vts/edges/faces
    t_get_inside.start();
    vector<int> vts_indices, edges_indices, faces_indices;
    getValidVts(vts_indices);
    getValidEdges(edges_indices);
    getValidFaces(faces_indices);
    t_get_inside.stop();
    // std::cout << "inside indices ready." << std::endl;

    // obtain the requested measure
    t_msure.start();
    std::cout << "Computing measures ..." << std::endl;
    computeVertexMeasure(MeasureForMA::LAMBDA, vts_indices, _vts_msure);
    computeEdgesMeasure(MeasureForMA::LAMBDA, edges_indices, _edges_msure);
    computeFacesMeasure(MeasureForMA::LAMBDA, faces_indices, _faces_msure);
    std::cout << "Done: computing measures" << std::endl;
    t_msure.stop();
    /*
    t_msure.start();
    _vts_msure.clear();
    for ( auto i : vts_indices )
            _vts_msure.push_back( m_v_msure[ i ] );
    _edges_msure.clear();
    for ( auto i : edges_indices )
            _edges_msure.push_back( m_e_msure[ i ] );
    _faces_msure.clear();
    for ( auto i : faces_indices )
            _faces_msure.push_back( m_f_msure[ i ] );
    t_msure.stop();
    */
    std::cout << "measures extracted." << std::endl;
    auto v_msure_range =
        std::minmax_element(_vts_msure.begin(), _vts_msure.end());
    auto e_msure_range =
        std::minmax_element(_edges_msure.begin(), _edges_msure.end());
    auto f_msure_range =
        std::minmax_element(_faces_msure.begin(), _faces_msure.end());
    std::cout << "V measure: [" << *v_msure_range.first << ","
              << *v_msure_range.second << "]" << std::endl;
    std::cout << "E measure: [" << *e_msure_range.first << ","
              << *e_msure_range.second << "]" << std::endl;
    std::cout << "F measure: [" << *f_msure_range.first << ","
              << *f_msure_range.second << "]" << std::endl;

    t_get_inside.restart();
    // grab the geometry of those elements
    getVts(vts_indices, _output_vts);
    // grab edges and keep the order (to merge in new edges later)
    getEdges(edges_indices, _output_edges);
    std::unordered_set<ivec2, ivec2Hash> unique_edge_set;
    for (auto ei = 0; ei < _output_edges.size(); ++ei)
        unique_edge_set.insert(_output_edges[ei]);
    vector<int> f_vrep;
    vector<uTriFace> f_tris;
    vector<float> f_msure_temp;
    ivec2 tri_edges[3];
    for (auto fit = faces_indices.begin(); fit != faces_indices.end(); fit++)
    {
        geom().getFaceVRep(*fit, f_vrep);
        util::simpleTriangulate(f_vrep, f_tris);
        float msure = _faces_msure[fit - faces_indices.begin()];
        for (auto i = 0; i < f_tris.size(); ++i)
        {
            f_msure_temp.push_back(msure);
            _output_tris.push_back(f_tris[i]);
            _from_fi.push_back(*fit);
            util::makeEdgesFromTri(f_tris[i], tri_edges);
            // add new edges arising from face being triangulated
            for (auto e : tri_edges)
            {
                if (unique_edge_set.count(e) == 0)
                {
                    _output_edges.push_back(e);
                    _edges_msure.push_back(msure); // same as the face's measure
                    unique_edge_set.insert(e);
                }
            }
        }
    }
    _faces_msure.swap(f_msure_temp);
    {
        auto tmp = decltype(f_msure_temp)();
        f_msure_temp.swap(tmp);
    }
    unique_edge_set.clear();
    t_get_inside.stop();

    std::cout << "inside geometry obtained." << std::endl;
    t_get_inside.restart();
    util::compactify(vts_indices, _output_edges, _output_tris);
    t_get_inside.stop();
    std::cout << "inside geometry compact-ified." << std::endl;

    cout << "time -> extract inside: " << t_get_inside.elapseMilli().count()
         << " ms" << endl;
    cout << "time -> compute measure: " << t_msure.elapseMilli().count()
         << " ms" << endl;
}

// void VoroInfo::mergeCloseVts( float _eps, bool _only_inside )
//{
//	/*timers to break down times into each step*/
//	timer temp_timer;
//	auto build_knn_t = timer::duration::zero();
//	auto search_knn_t = timer::duration::zero();
//	auto merge_opt_t = timer::duration::zero();
//	auto finalize_t = timer::duration::zero();
//	/* maintain important info */
//	// for vertex: before merging index -> after merging index (into the set
// of remaining vts) 	vector<int> old_to_new_idx;
//	// the remaining vertices after merging (indices into the old vts list)
//	vector<int> new_vts_indices;
//	// the set of remaining edges which might be important to preserving
// topology after merging 	unordered_set<ivec2, ivec2Hash> new_edges;
//	// remove tag for each face
//	vector<bool> remove_face;
//	// the cell complex representing new geometry after merging
//	cellcomplex new_geom;
//	new_geom.allocVts( m_geom.numVts() / 2 );
//	new_geom.allocEdges( m_geom.numEdges() / 2 );
//	new_geom.allocFaces( m_geom.numFaces() / 2 );
//	// debug
//	//set<int> old_vts_dbg = { 332, 294 };

//	temp_timer.start();
//	float s = 1000000.0f;
//	// insert voro vts into kdtree
//	const auto& vts = m_geom.getAllVts();
//	int n_pts = vts.size();
//	ANNpointArray data_pts;
//	ANNpoint query_p;
//	ANNidxArray nn_idx;
//	ANNkd_tree* kdtree;
//	query_p = annAllocPt( 3 );
//	data_pts = annAllocPts( n_pts, 3 );
//	for ( size_t i = 0; i < n_pts; ++i )
//	{
//		/*if ( old_vts_dbg.count( i ) )
//			int stop = 1;*/
//		auto v = s * vts[ i ];
//		data_pts[ i ][ 0 ] = v[ 0 ];
//		data_pts[ i ][ 1 ] = v[ 1 ];
//		data_pts[ i ][ 2 ] = v[ 2 ];
//	}
//	kdtree = new ANNkd_tree( data_pts, n_pts, 3 );
//	temp_timer.stop();
//	build_knn_t = temp_timer.elapse();
//	cout << "done: inserting vertex data to kdtree." << endl;

//	// Merge vertices:
//	// go over each vert, find all close vts, "merge" them to one vertex
//	vector<bool> visited;
//	visited.assign( n_pts, false );
//	old_to_new_idx.resize( n_pts, -1 );
//	float sq_r = _eps * _eps; cout << "sq. radius = " << sq_r << endl;
//	int nn_within_r = 1;
//	for ( size_t vi = 0; vi < n_pts; ++vi )
//	{
//		/*if ( old_vts_dbg.count( vi ) )
//			int stop = 1;*/
//		if ( visited[ vi ] )
//			continue;

//		if ( _only_inside && !isVertexValid(vi) )
//			continue;

//		new_vts_indices.push_back( vi );
//		auto v = s * vts[ vi ];
//		query_p[ 0 ] = v[ 0 ];
//		query_p[ 1 ] = v[ 1 ];
//		query_p[ 2 ] = v[ 2 ];

//		// obtain all nearest neighbors within the ball of the specified
// radius 		temp_timer.start(); 		nn_within_r =
// kdtree->annkFRSearch( query_p, sq_r, 0, nullptr, nullptr, 0.0 );
// if ( nn_within_r == 0 ) 			cout << "WARNING: vert
//"<<vi << ", zero NN found!" << endl; 		nn_idx = new ANNidx[ nn_within_r
//]; 		kdtree->annkFRSearch( query_p, sq_r, nn_within_r, nn_idx,
// nullptr, 0.0 );
// temp_timer.stop(); 		search_knn_t += temp_timer.elapse();
//
//		// assign to those n.n. new index
//		temp_timer.start();
//		for ( size_t j = 0; j < nn_within_r; ++j )
//		{
//			auto idx = nn_idx[ j ];
//			visited[ idx ] = true;
//			old_to_new_idx[ idx ] = new_vts_indices.size() - 1;
//		}
//		temp_timer.stop();
//		merge_opt_t += temp_timer.elapse();
//		delete [] nn_idx;
//	}
//	cout << "done: merging close vertices." << endl;
//	//for ( size_t i = 0; i < old_to_new_idx.size(); ++i )
//	//{
//	//	//cout << "old " << i << " -> new " << old_to_new_idx[ i ] <<
// endl;
//	//	cout << old_to_new_idx[ i ] << " ";
//	//}
//	//cout << endl;
//	temp_timer.start();
//	delete kdtree;
//	annClose();
//	temp_timer.stop();
//	build_knn_t += temp_timer.elapse();

//	temp_timer.start();
//	// Save remaining vts to the geometry struct
//	for ( const auto& vi : new_vts_indices )
//	{
//		new_geom.appendVert( m_geom.getVert( vi ) );
//	}
//	cout << "done: saving remaining vertices (after merging) into geometry
// struct." << endl;

//	// update validity for each new vertex
//	vector<bool> new_v_valid_tag;
//	vector<bool> new_v_finite_tag;
//	for ( const auto& vi : new_vts_indices )
//	{
//		new_v_valid_tag.push_back( m_vts_valid[ vi ] );
//		new_v_finite_tag.push_back( isVertexFinite( vi ) );
//	}
//	m_vts_valid = new_v_valid_tag;
//	m_is_finite_v = new_v_finite_tag;
//	new_v_valid_tag.clear();
//	new_v_valid_tag.shrink_to_fit();
//	new_v_finite_tag.clear();
//	new_v_finite_tag.shrink_to_fit();
//	cout << "done: updating validity tag for vertices." << endl;
//	// update radius function for each new vertex
//	if ( radiiValid() )
//	{
//		vector<float> new_r_per_v;
//		for ( const auto& vi : new_vts_indices )
//		{
//			new_r_per_v.push_back( m_r_per_v[ vi ] );
//		}
//		m_r_per_v = new_r_per_v;
//		new_r_per_v.clear();
//		new_r_per_v.shrink_to_fit();
//		cout << "done: updating radius function for vertices." << endl;
//	}

//	// Checking edges:
//	// some edges may degenerate to vertex, some may still are valid edges.
//	for ( size_t ei = 0; ei < m_geom.numEdges(); ++ei )
//	{
//		if ( _only_inside && !isEdgeValid( ei ) )
//			continue;

//		const auto& e = m_geom.getEdge( ei );
//		//cout << "old edge " << ei << ": " << e << endl;
//		auto u = old_to_new_idx[ e[ 0 ] ];
//		auto v = old_to_new_idx[ e[ 1 ] ];
//		auto new_e = util::makeEdge( u, v );

//		if ( u == v )// invalid edge. skip.
//			continue;
//
//		//cout << "new edge " << new_e << endl;
//		new_edges.insert( new_e );
//	}
//	cout << "done: checking and handling edges degeneracy after merging." <<
// endl;

//	// Checking faces:
//	// some faces may degenerate to edges, some may still be valid faces.
//	// update them accordingly. associated attributes, like contributing
// sites
//	// to those faces may need to be merged to near-by faces' attributes
// list 	auto n_faces = m_geom.numFaces(); 	remove_face.assign(
// n_faces,
// false ); 	unordered_set<int> unique_vts_of_new_f; 	vector<int>
// vts_of_old_f;
//	vector<int> vts_of_new_f;
//	int cur_v;
//	for ( size_t fi = 0; fi < n_faces; ++fi )
//	{
//		/*if ( fi == 244 )
//			int stop = 1;*/
//		if ( _only_inside && !isFaceValid( fi ) )
//		{
//			remove_face[ fi ] = true;
//			continue;
//		}

//		m_geom.getFaceVRep( fi, vts_of_old_f );
//		unique_vts_of_new_f.clear();
//		vts_of_new_f.clear();
//		for ( size_t j = 0; j < vts_of_old_f.size(); ++j )
//		{
//			cur_v = old_to_new_idx[ vts_of_old_f[ j ] ];
//			if (unique_vts_of_new_f.count(cur_v) > 0 )
//				continue;
//			unique_vts_of_new_f.insert( cur_v );
//			vts_of_new_f.push_back( cur_v );
//		}

//		//// reject if this is not an *inside face*
//		//if ( !isFaceValid( vts_of_new_f ) )
//		//	continue;

//		// check how degenerate this face has become
//		if ( vts_of_new_f.size() < 2 )
//		{ // degenerate to a vertex. remove it.
//			remove_face[ fi ] = true;
//		}
//		else if ( vts_of_new_f.size() == 2 )
//		{ // degenerate to an edge. add to edge-list
//			new_edges.insert( util::makeEdge( vts_of_new_f[ 0 ],
// vts_of_new_f[ 1 ] ) ); 			remove_face[ fi ] = true;
//		}
//		else // ( vts_of_new_f.size() > 2 )
//		{ // still a valid face. add to face-list
//			new_geom.appendFace( vts_of_new_f );
//			remove_face[ fi ] = false;

//			// debug
//			/*if ( vts_of_new_f.size() < vts_of_old_f.size() )
//			{
//				cout << "face " << fi << ": " << endl;
//				cout << "old ";
//				for ( auto j : vts_of_old_f )
//					cout << j << m_geom.getVert( j ) << ",";
//				cout << endl;
//				cout << "new ";
//				for ( auto j : vts_of_new_f )
//					cout << j << new_geom.getVert( j ) <<
//","; 				cout << endl;
//			}*/
//		}
//		/*if ( !vts_of_new_f.empty() )
//		{
//			cout << "old: ";
//			for ( auto v : vts_of_old_f )
//			{
//				cout << v << " ";
//			}
//			cout << ",";
//			cout << "new: ";
//			for ( auto v : vts_of_new_f )
//			{
//				cout << v << " ";
//			}
//			cout << ",";
//		}*/
//	}
//	//cout << endl;
//	cout << "done: checking and handling faces degeneracy after merging." <<
// endl;

//	// Save new edges (important to topo preserving) to geometry struct
//	for ( const auto& e : new_edges )
//	{
//		new_geom.appendEdge( e );
//	}
//	new_edges.clear();
//	cout << "done: saving edges into geometry struct." << endl;
//	temp_timer.stop();
//	merge_opt_t += temp_timer.elapse();

//	/*Finalize the geometry struct*/
//	temp_timer.start();
//	new_geom.finalize();
//	new_geom.collectSpace();
//	/*for ( size_t ei = 0; ei < new_geom.numEdges(); ++ei )
//	{
//		if ( new_geom.refCount( ei ) == 0 )
//			cout << "edge " << ei << " refcnt == 0" << endl;
//	}*/

//	// replace old geometry with the new one
//	m_geom = new_geom;
//	new_geom.clear();
//	cout << "done: replacing old with new geometry struct." << endl;

//	/* update edge and face flags */
//	this->recomputeEdgeFlags();
//	this->recomputeFaceFlags();

//	/* update voro-face-related attributes */
//	// update sites:
//	// here we remove all sites for those faces are marked true in
// remove_face[] 	vector<ivec2> facesites_cpy; 	for ( size_t i = 0; i <
// remove_face.size(); ++i )
//	{
//		if ( !remove_face[ i ] )
//			facesites_cpy.push_back( m_face_sites[ i ] );
//	}
//	m_face_sites = facesites_cpy;
//	facesites_cpy.clear(); facesites_cpy.shrink_to_fit();
//	// update edge / face is-finite-flag
//
//	temp_timer.stop();
//	finalize_t = temp_timer.elapse();
//	cout << "done: updating sites for voro-faces." << endl;

//	cout << "-------times break-down-------" << endl;
//	using namespace std::chrono;
//	cout << "knn tree build & destroy overhead: "
//		<< duration_cast<milliseconds>( build_knn_t ).count() << "ms" <<
// endl; 	cout << "knn search: "
//		<< duration_cast<milliseconds>(search_knn_t).count() << "ms" <<
// endl; 	cout << "merging operation: "
//		<< duration_cast<milliseconds>(merge_opt_t).count() << "ms" <<
// endl; 	cout << "voro finalize: "
//		<< duration_cast<milliseconds>(finalize_t).count() << "ms" <<
// endl; 	cout << "-------(times break-down)-------" << endl;
//}

bool VoroInfo::generateMeasure(MeasureForMA::meassuretype _msure_type)
{
    if (!geom().isFinalized())
        return false;
    m_v_msure.clear();
    m_e_msure.clear();
    m_f_msure.clear();
    vector<int> all_V(geom().numVts());
    std::iota(all_V.begin(), all_V.end(), 0);
    computeVertexMeasure(MeasureForMA::LAMBDA, all_V, m_v_msure);
    all_V.resize(geom().numEdges());
    std::iota(all_V.begin(), all_V.end(), 0);
    computeEdgesMeasure(MeasureForMA::LAMBDA, all_V, m_e_msure);
    all_V.resize(geom().numFaces());
    std::iota(all_V.begin(), all_V.end(), 0);
    computeFacesMeasure(MeasureForMA::LAMBDA, all_V, m_f_msure);

    return true;
}

void VoroInfo::mergeCloseVts(float _eps, bool _only_inside)
{
    /*timers to break down times into each step*/
    timer temp_timer;
    auto merge_opt_t = timer::duration::zero();
    auto finalize_t = timer::duration::zero();
    /* maintain important info */
    // for vertex: before merging index -> after merging index (into the set of
    // remaining vts)
    vector<int> old_to_new_idx;
    // the remaining vertices after merging (indices into the old vts list)
    vector<int> new_vts_indices;
    // the set of remaining edges which might be important to preserving
    // topology after merging
    unordered_set<ivec2, ivec2Hash> new_edges;
    // remove tag for each face
    vector<bool> remove_face;
    // the cell complex representing new geometry after merging
    cellcomplex new_geom;
    new_geom.allocVts(m_geom.numVts() / 2);
    new_geom.allocEdges(m_geom.numEdges() / 2);
    new_geom.allocFaces(m_geom.numFaces() / 2);
    // debug
    // set<int> old_vts_dbg = { 332, 294 };

    int n_pts = geom().numVts();

    /* Merge vertices: */
    temp_timer.start();
    vector<bool> visited;
    queue<int> close_q;
    visited.assign(n_pts, false);
    old_to_new_idx.resize(n_pts, -1);
    float sq_r = _eps * _eps;
    cout << "sq. radius = " << sq_r << endl;
    int nn_within_r = 1;
    // go over each vert, find all close vts, "merge" them to one vertex
    for (size_t vi = 0; vi < n_pts; ++vi)
    {
        /*if ( old_vts_dbg.count( vi ) )
        int stop = 1;*/
        if (visited[vi])
            continue;

        if (_only_inside && !isVertexValid(vi))
            continue;

        new_vts_indices.push_back(vi);

        // merge all nbs that are close to vi
        geom().needVVAdjacency();
        close_q.push(vi);
        visited[vi] = true;
        old_to_new_idx[vi] = new_vts_indices.size() - 1;
        vector<int> nb_vs;
        nb_vs.reserve(50);
        while (!close_q.empty())
        {
            auto cur = close_q.front();
            close_q.pop();
            nb_vs.clear();
            geom().nbVts(cur, nb_vs);
            for (auto ni : nb_vs)
            {
                if (!visited[ni] && trimesh::dist2(geom().getVert(vi),
                                                   geom().getVert(ni)) <= sq_r)
                {
                    close_q.push(ni);
                    visited[ni] = true;
                    old_to_new_idx[ni] = new_vts_indices.size() - 1;
                }
            }
        }
    }
    temp_timer.stop();
    merge_opt_t += temp_timer.elapse();
    cout << "done: merging close vertices." << endl;
    // for ( size_t i = 0; i < old_to_new_idx.size(); ++i )
    //{
    //	//cout << "old " << i << " -> new " << old_to_new_idx[ i ] << endl;
    //	cout << old_to_new_idx[ i ] << " ";
    // }
    // cout << endl;

    temp_timer.start();
    // Save remaining vts to the geometry struct
    for (const auto& vi : new_vts_indices)
    {
        new_geom.appendVert(m_geom.getVert(vi));
    }
    cout << "done: saving remaining vertices (after merging) into geometry "
            "struct."
         << endl;

    // update validity for each new vertex
    if (!onlyHasInside())
    {
        vector<bool> new_v_valid_tag;
        vector<bool> new_v_finite_tag;
        for (const auto& vi : new_vts_indices)
        {
            new_v_valid_tag.push_back(m_vts_valid[vi]);
            new_v_finite_tag.push_back(isVertexFinite(vi));
        }
        m_vts_valid = new_v_valid_tag;
        m_is_finite_v = new_v_finite_tag;
        new_v_valid_tag.clear();
        new_v_valid_tag.shrink_to_fit();
        new_v_finite_tag.clear();
        new_v_finite_tag.shrink_to_fit();
        cout << "done: updating validity tag for vertices." << endl;
    }
    // update radius function for each new vertex
    if (radiiValid())
    {
        vector<float> new_r_per_v;
        for (const auto& vi : new_vts_indices)
        {
            new_r_per_v.push_back(m_r_per_v[vi]);
        }
        m_r_per_v = new_r_per_v;
        new_r_per_v.clear();
        new_r_per_v.shrink_to_fit();
        cout << "done: updating radius function for vertices." << endl;
    }
    // keep only those vert-msure for the new vts
    if (isMsureReady())
    {
        vector<float> new_v_msure;
        for (auto i : new_vts_indices)
            new_v_msure.push_back(m_v_msure[i]);
        m_v_msure = std::move(new_v_msure);
        cout << "done: updating measure for vertices." << endl;
    }

    // Checking edges:
    // some edges may degenerate to vertex, some may still are valid edges.
    // We also want to keep only the measures for the valid edges.
    // NOTE: We want to make sure the correspondence between each edge and its
    // measure is correct This could be a concern because later the new-geom is
    // going to be "finalized" which will rebuild edges list (from the set of
    // edges used by all faces) HOWEVER, the rebuilt edge list SHOULD be the
    // same as the new edges extracted here
    vector<float> new_e_msure; // (keep those valid edges' measure)
    for (size_t ei = 0; ei < m_geom.numEdges(); ++ei)
    {
        if (_only_inside && !isEdgeValid(ei))
            continue;

        const auto& e = m_geom.getEdge(ei);
        // cout << "old edge " << ei << ": " << e << endl;
        auto u = old_to_new_idx[e[0]];
        auto v = old_to_new_idx[e[1]];

        if (u == v) // invalid edge. skip.
            continue;

        auto new_e = util::makeEdge(u, v);
        // cout << "new edge " << new_e << endl;
        new_edges.insert(new_e);
        if (isMsureReady())
            new_e_msure.push_back(m_e_msure[ei]); // save msure for this edge
    }
    if (isMsureReady())
        m_e_msure = std::move(new_e_msure);
    cout << "done: checking and handling edges degeneracy after merging (also "
            "might have kept valid edges' measure)."
         << endl;

    // Checking faces:
    // some faces may degenerate to edges, some may still be valid faces.
    // update them accordingly. associated attributes, like contributing sites
    // to those faces may need to be merged to near-by faces / edges
    auto n_faces = m_geom.numFaces();
    remove_face.assign(n_faces, false);
    unordered_set<int> unique_vts_of_new_f;
    vector<int> vts_of_old_f;
    vector<int> vts_of_new_f;
    int cur_v;
    for (size_t fi = 0; fi < n_faces; ++fi)
    {
        /*if ( fi == 244 )
        int stop = 1;*/
        if (_only_inside && !isFaceValid(fi))
        {
            remove_face[fi] = true;
            continue;
        }

        m_geom.getFaceVRep(fi, vts_of_old_f);
        unique_vts_of_new_f.clear();
        vts_of_new_f.clear();
        for (size_t j = 0; j < vts_of_old_f.size(); ++j)
        {
            cur_v = old_to_new_idx[vts_of_old_f[j]];
            if (unique_vts_of_new_f.count(cur_v) > 0)
                continue;
            unique_vts_of_new_f.insert(cur_v);
            vts_of_new_f.push_back(cur_v);
        }

        //// reject if this is not an *inside face*
        // if ( !isFaceValid( vts_of_new_f ) )
        //	continue;

        // check how degenerate this face has become
        if (vts_of_new_f.size() < 2)
        { // degenerate to a vertex. remove it.
            remove_face[fi] = true;
        }
        else if (vts_of_new_f.size() == 2)
        { // degenerate to an edge. add to edge-list
            new_edges.insert(util::makeEdge(vts_of_new_f[0], vts_of_new_f[1]));
            remove_face[fi] = true;
        }
        else // ( vts_of_new_f.size() > 2 )
        {    // still a valid face. add to face-list
            new_geom.appendFace(vts_of_new_f);
            remove_face[fi] = false;

            // debug
            /*if ( vts_of_new_f.size() < vts_of_old_f.size() )
            {
            cout << "face " << fi << ": " << endl;
            cout << "old ";
            for ( auto j : vts_of_old_f )
            cout << j << m_geom.getVert( j ) << ",";
            cout << endl;
            cout << "new ";
            for ( auto j : vts_of_new_f )
            cout << j << new_geom.getVert( j ) << ",";
            cout << endl;
            }*/
        }
        /*if ( !vts_of_new_f.empty() )
        {
        cout << "old: ";
        for ( auto v : vts_of_old_f )
        {
        cout << v << " ";
        }
        cout << ",";
        cout << "new: ";
        for ( auto v : vts_of_new_f )
        {
        cout << v << " ";
        }
        cout << ",";
        }*/
    }
    // cout << endl;
    cout << "done: checking and handling faces degeneracy after merging."
         << endl;

    // Save new edges (important to topo preserving) to geometry struct
    for (const auto& e : new_edges)
    {
        new_geom.appendEdge(e);
    }
    new_edges.clear();
    cout << "done: saving edges into geometry struct." << endl;
    temp_timer.stop();
    merge_opt_t += temp_timer.elapse();

    /*Finalize the new geometry struct*/
    temp_timer.start();
    //// estimate circumradius for edges w/o nb valid faces
    //// (this is hacky but efficient. A more consistent but slower way is we
    /// should associate the sites / for the nb invalid faces to the edge, then
    /// compute circum-radius later when lambda measure is requested)
    // for ( size_t ei = 0; ei < new_geom.numEdges(); ++ei )
    //{
    //	if ( new_geom.refCount( ei ) == 0 )
    //	{

    //	}
    //}

    // replace old geometry with the new one
    m_geom = std::move(new_geom);
    m_geom.finalize();
    m_geom.collectSpace();
    // new_geom.clear();
    cout << "done: replacing old with new geometry struct." << endl;

    /* update edge and face flags */
    this->recomputeEdgeFlags();
    this->recomputeFaceFlags();

    /* update voro-face-related attributes */
    // per-face sites & measure
    // here we remove all sites for those faces are marked true in remove_face[]
    vector<ivec2> facesites_cpy;
    vector<float> new_f_msure;
    for (size_t i = 0; i < remove_face.size(); ++i)
    {
        if (!remove_face[i])
        {
            facesites_cpy.push_back(m_face_sites[i]);
            if (isMsureReady())
                new_f_msure.push_back(m_f_msure[i]);
        }
    }
    m_face_sites = std::move(facesites_cpy);
    if (isMsureReady())
        m_f_msure = std::move(new_f_msure);
    // update edge / face is-finite-flag

    temp_timer.stop();
    finalize_t = temp_timer.elapse();
    cout << "done: updating sites & measure for faces." << endl;

    // cout << "-------times break-down-------" << endl;
    using namespace std::chrono;
    cout << "time -> merging operation: "
         << duration_cast<milliseconds>(merge_opt_t).count() << " ms" << endl;
    /*cout << "time -> voro finalize: "
            << duration_cast<milliseconds>( finalize_t ).count() << " ms" <<
       endl;*/
    // cout << "-------(times break-down)-------" << endl;
}

void VoroInfo::computeEulerChar(eulerchar& _euler_struct, bool _do_prune_first)
{
    // compute_multiplicity_vts();
    int V = 0, E = 0, F = 0, C = 0, euler;
    int V_s = 0;
    vector<int> f_v_rep;
    vector<int> f_e_rep;
    if (_do_prune_first)
    {
        // vts remained after pruning
        vector<int> valid_vts_indices;
        getValidVts(valid_vts_indices);
        // edges remained after pruning
        vector<int> valid_edges_indices;
        getValidEdges(valid_edges_indices);
        // store faces lying inside the shape
        vector<int> valid_faces_indices;
        getValidFaces(valid_faces_indices);
        V = valid_vts_indices.size();
        E = (int)valid_edges_indices.size();
        F = (int)valid_faces_indices.size();
        // TODO: need to take in a set of valid vts, not faces
        C = m_geom.compute_conn_cmpnts(&valid_faces_indices, 1);
    }
    else
    {
        V = (int)m_geom.numVts();
        E = (int)m_geom.numEdges();
        F = (int)m_geom.numFaces();
        C = m_geom.compute_conn_cmpnts();
    }

    // for ( const auto& vi : remain_vts_indices )
    //{
    //	auto mult_v = m_mult_num_vts[ vi ];
    //	if ( mult_v > 1 )
    //		V_s++;
    //	//cout << "mult of v: " << mult_v << endl;
    // }
    // std::cout << "# singular V = " << V_s << std::endl;

    euler = V - E + F;

    _euler_struct = eulerchar(V, E, F, 0, C, euler);
}

void VoroInfo::computeVertexMeasure(MeasureForMA::meassuretype _mssure_tp,
                                    const vector<int>& _vts_indices,
                                    vector<float>& _vts_msure) const
{
    _vts_msure.clear();
    _vts_msure.reserve(_vts_msure.size());
    float v_lmd = 0.0f;
    switch (_mssure_tp)
    {
            // lambda of a vertex should be the circumradius
            // of the sites associated with the nb faces, which is
            // given by the radius at the vertex
        case MeasureForMA::LAMBDA:
        {
            vector<point> foot_pts;
            for (auto it = _vts_indices.begin(); it != _vts_indices.end(); ++it)
            {
                // collect all nearby sites
                foot_pts.clear();
                v_lmd = 0.0f;
                int n_nb_e = geom().cntNbEdgesofVert(*it);
                for (auto ni = 0; ni < n_nb_e; ++ni)
                {
                    auto nb_ei = geom().nbEdgeofVert(*it, ni);
                    int n_nb_f = geom().cntNbFacesofEdge(nb_ei);
                    for (auto nj = 0; nj < n_nb_f; ++nj)
                    {
                        auto nb_fi = geom().nbFaceofEdge(nb_ei, nj);
                        if (!isFaceValid(nb_fi))
                            continue;
                        const auto& sites = getSitesOfFace(nb_fi);
                        auto f_lmd = MeasureForMA::lambdaForFace(
                            m_site_positions[sites[0]],
                            m_site_positions[sites[1]]);
                        v_lmd = std::max(v_lmd, f_lmd);
                        /*foot_pts.push_back( m_site_positions[ sites[ 0 ] ] );
                        foot_pts.push_back( m_site_positions[ sites[ 1 ] ] );*/
                    }
                }
                /*if ( foot_pts.empty() )
                        std::cout << "WARNING: 0 nearby sites obtained for voro
                   v " << *it
                   << std::endl;*/
                // compute an approximate circumradius of all those sites
                /*point cent = trimesh::point_center_of_mass( foot_pts );
                v_lmd = 0.0f;
                for ( const auto& p : foot_pts )
                        v_lmd = std::max( v_lmd, trimesh::dist( cent, p ) );*/
                _vts_msure.push_back(v_lmd);
            }
            break;
        }
        default:
            v_lmd = 0.0f;
            break;
    }
}

void VoroInfo::computeEdgesMeasure(MeasureForMA::meassuretype _mssure_tp,
                                   const vector<int>& _edges_indices,
                                   vector<float>& _edges_msure) const
{
    _edges_msure.clear();
    _edges_msure.reserve(_edges_indices.size());
    switch (_mssure_tp)
    {
        case MeasureForMA::LAMBDA:
        {
            vector<point> foot_pts;
            float e_lmd = 0.0f;
            for (auto it = _edges_indices.begin(); it != _edges_indices.end();
                 it++)
            {
                e_lmd = 0.0f;
                foot_pts.clear();
                int n_nb_f = geom().cntNbFacesofEdge(*it);
                if (n_nb_f > 0)
                {
                    for (auto j = 0; j < n_nb_f; ++j)
                    {
                        auto nb_fi = geom().nbFaceofEdge(*it, j);
                        if (!isFaceValid(nb_fi))
                            continue;
                        const auto& sites = getSitesOfFace(nb_fi);
                        auto f_lmd = MeasureForMA::lambdaForFace(
                            m_site_positions[sites[0]],
                            m_site_positions[sites[1]]);
                        e_lmd = std::max(e_lmd, f_lmd);
                        /*foot_pts.push_back( m_site_positions[ sites[ 0 ] ] );
                        foot_pts.push_back( m_site_positions[ sites[ 1 ] ] );*/
                    }
                }
                /*if ( foot_pts.empty() )
                        std::cout << "WARNING: 0 nearby sites obtained for voro
                edge " << *it << std::endl; auto cent =
                trimesh::point_center_of_mass( foot_pts );*/
                /*e_lmd = 0.0f;
                for ( const auto& p : foot_pts )
                        e_lmd = std::max( e_lmd, trimesh::dist( cent, p ) );*/
                _edges_msure.push_back(e_lmd);
            }
            break;
        }
        default:
            break;
    }
}
void VoroInfo::computeEdgesMeasure(const vector<float>& _vts_msure,
                                   const vector<int>& _edges_indices,
                                   vector<float>& _edges_msure) const
{
    _edges_msure.clear();
    _edges_msure.reserve(_edges_indices.size());
    for (auto i = 0; i < _edges_indices.size(); ++i)
    {
        auto e = geom().getEdge(_edges_indices[i]);
        _edges_msure.push_back((_vts_msure[e[0]] + _vts_msure[e[1]]) * 0.5f);
    }
}

void VoroInfo::computeFacesMeasure(MeasureForMA::meassuretype _mssure_tp,
                                   const vector<int>& _faces_indices,
                                   vector<float>& _faces_msure) const
{
    _faces_msure.clear();
    _faces_msure.reserve(_faces_indices.size());
    switch (_mssure_tp)
    {
        case MeasureForMA::LAMBDA:
            for (auto it = _faces_indices.begin(); it != _faces_indices.end();
                 it++)
            {
                // lambda of a face is the distance between two sites
                const auto& sites = getSitesOfFace(*it);
                float f_lambda = MeasureForMA::lambdaForFace(
                    m_site_positions[sites[0]], m_site_positions[sites[1]]);
                _faces_msure.push_back(f_lambda);
            }
            break;
        default:
            break;
    }
}

void VoroInfo::outputToMathematica(const char* _filename) const
{
    ofstream ofile(_filename);
    if (!ofile.is_open())
    {
        cout << "Error: cannot open file: " << _filename << endl;
        return;
    }

    cout << "Writing current voro geometry (2-complex) to Mathematica..."
         << endl;

    // output outmost {
    ofile << '{' << endl;

    // output vertices (all)
    ofile << "{(* begin: vertices *)" << endl;
    auto default_prec = ofile.precision();
    ofile << std::fixed << std::setprecision(8);
    for (size_t i = 0; i < m_geom.numVts(); ++i)
    {
        auto v = m_geom.getVert(i);
        ofile << v[0] << ',' << v[1] << ',' << v[2];
        if (i < m_geom.numVts() - 1)
        {
            ofile << ",";
        }
    }
    ofile << "},"
          << "(* end: vertices *)" << endl;
    ofile << std::defaultfloat << std::setprecision(default_prec);

    // output edges (0-ref.ed, finite)
    ofile << "{(* begin: edges *)" << endl;
    for (size_t i = 0; i < m_geom.numEdges(); ++i)
    {
        if (!isEdgeValid(i) || m_geom.refCount(i) > 0)
            continue;
        auto e = m_geom.getEdge(i);
        ofile << e[0] + 1 << "," << e[1] + 1 << ",";
    }
    ofile.seekp(-1, ofile.cur);
    ofile << "},"
          << "(* end: edges *)" << endl;

    // output faces (finite)
    ofile << "{(* begin: faces *)" << endl;
    vector<int> vts_of_f;
    for (size_t i = 0; i < m_geom.numFaces(); ++i)
    {
        if (!isFaceValid(i))
            continue;
        vts_of_f.clear();
        m_geom.getFaceVRep(i, vts_of_f);
        ofile << "{";
        for (int j = 0; j < vts_of_f.size(); ++j)
        {
            ofile << vts_of_f[j] + 1;
            if (j < vts_of_f.size() - 1)
                ofile << ",";
        }
        ofile << "},";
    }
    ofile.seekp(-1, ofile.cur);
    ofile << "}"
          << "(* end: faces *)" << endl;

    // output outmost }
    ofile << '}' << endl;
    ofile.close();

    cout << "Done: Writing current voro geometry (2-complex) to Mathematica."
         << endl;
}

void VoroInfo::invalidate_geometric_states()
{
    m_r_valid = false;
    m_face_sites_valid = false;
    m_only_inside = -1;
}

void VoroInfo::set_only_inside(bool yesno) { m_only_inside = yesno; }

void VoroInfo::load_voro_vts(ifstream& _in_node,
                             shared_ptr<Volume3DScalar> _vol,
                             unordered_map<int, int>* _old_new_v_map)
{
    // read in content all at once and parse numbers using a stringstream
    _in_node.seekg(0, std::ios::end);
    auto len = _in_node.tellg();
    _in_node.seekg(std::ios::beg);
    vector<char> buffer(len);
    _in_node.read(&buffer[0], len);

    m_geom.clearVertList();
    this->m_is_finite_v.clear();
    int n_vts;
    int nchar_read = 0, offset = 0;
    auto delim = " \n\t";

    auto token = strtok(&buffer[0], delim);
    n_vts = atoi(token);
    // skip 3 dump ints
    token = strtok(NULL, delim);
    token = strtok(NULL, delim);
    token = strtok(NULL, delim);

    point v;
    int vid;
    m_bbox.clear();
    for (size_t i = 0; i < n_vts; ++i)
    {
        // skip the node id
        token = strtok(NULL, delim);
        vid = atoi(token);
        assert(vid == i);
        // read x, y, z of node's coords
        token = strtok(NULL, delim);
        v[0] = atof(token);
        token = strtok(NULL, delim);
        v[1] = atof(token);
        token = strtok(NULL, delim);
        v[2] = atof(token);

        if (_vol && this->tagVert(v, _vol) == false)
            continue;
        if (_vol)
            (*_old_new_v_map)[i] = m_geom.numVts();
        m_geom.appendVert(v);
        if (!_vol) // only need *finite* when loading entire VD
            m_is_finite_v.push_back(true);
        m_bbox += v;
    }

    buffer.clear();
    buffer.shrink_to_fit();
}

void VoroInfo::load_voro_edges(
    ifstream& _in_edge, shared_ptr<Volume3DScalar> _vol,
    const unordered_map<int, int>* const _old_new_v_map,
    unordered_map<int, int>* _old_new_e_map)
{
    m_geom.clearEdgeList();
    m_is_finite_e.clear();

    _in_edge.seekg(0, std::ios::end);
    auto len = _in_edge.tellg();
    _in_edge.seekg(std::ios::beg);
    vector<char> buffer(len);
    _in_edge.read(&buffer[0], len);

    int n_edges, dump_int;
    auto delim = " \n\t";
    auto token = strtok(&buffer[0], delim);
    // num of edges info
    n_edges = atoi(token);
    // skip irrelevant int
    token = strtok(NULL, delim);

    // preallocate n_edges info
    m_is_finite_e.reserve(n_edges);

    for (int eid = 0; eid < n_edges; ++eid)
    {
        ivec2 e;
        trimesh::vec dir; // dir for infinite edge
        // skip edge id
        token = strtok(NULL, delim);
        // read edge info
        token = strtok(NULL, delim);
        e[0] = atoi(token);
        token = strtok(NULL, delim);
        e[1] = atoi(token);

        if (e[1] == -1) // edge is infinite.
        {
            // read its direction
            token = strtok(NULL, delim);
            dir[0] = atof(token);
            token = strtok(NULL, delim);
            dir[1] = atof(token);
            token = strtok(NULL, delim);
            dir[2] = atof(token);
        }

        if (_vol)
        {
            if (!_old_new_v_map->count(e[0]) || !_old_new_v_map->count(e[1]))
                continue;
            e[0] = _old_new_v_map->find(e[0])->second;
            e[1] = _old_new_v_map->find(e[1])->second;
            (*_old_new_e_map)[eid] = m_geom.numEdges();
        }

        // only handle element *finiteness* when loading entire VD
        if (!_vol)
            if (e[1] == -1)
            {
                m_is_finite_e.push_back(false);
                // create an infinite vertex for this infinite edge to refer to
                point inf_v =
                    geom().getVert(e[0]) + (-dir) * m_bbox.radius() * 0.01f;
                m_geom.appendVert(inf_v);
                m_is_finite_v.push_back(false);
                e[1] = m_geom.numVts() - 1;
            }
            else
            {
                m_is_finite_e.push_back(true);
            }
        m_geom.appendEdge(e);
    }
    m_is_finite_e.shrink_to_fit();

    buffer.clear();
    buffer.shrink_to_fit();
}

void VoroInfo::load_voro_faces(
    ifstream& _in_face, shared_ptr<Volume3DScalar> _vol,
    const unordered_map<int, int>* const _old_new_e_map)
{
    m_geom.clearFaceList();
    m_face_sites.clear();
    m_is_finite_f.clear();

    _in_face.seekg(0, std::ios::end);
    auto len = _in_face.tellg();
    _in_face.seekg(std::ios::beg);
    vector<char> buffer(len);
    _in_face.read(&buffer[0], len);

    int n_face;
    auto delim = " \n\t";
    // read in num of faces
    auto token = strtok(&buffer[0], delim);
    n_face = atoi(token);
    // skip irrelevant int
    token = strtok(NULL, delim);

    /*m_faceIndex.reserve( n_face );
    m_faceVtsN.reserve( n_face );*/
    m_face_sites.reserve(n_face);
    m_is_finite_f.reserve(n_face);

    vector<int> f_erep;
    int s1, s2, nsides, eid;
    vector<int> f_of_vts;
    for (size_t i = 0; i < n_face; ++i)
    {
        if (i == 505)
            int stop = 1;
        // skip face id
        token = strtok(NULL, delim);
        // read the two sites for this face
        token = strtok(NULL, delim);
        s1 = atoi(token);
        token = strtok(NULL, delim);
        s2 = atoi(token);
        // how many sides does this face have?
        token = strtok(NULL, delim);
        nsides = atoi(token);

        // read in nsides edges indices for this face
        f_erep.clear();
        for (size_t j = 0; j < nsides; ++j)
        {
            token = strtok(NULL, delim);
            eid = atoi(token);
            f_erep.push_back(eid);
        }

        // if vol is given, use it to do prune test
        if (_vol)
            for (size_t j = 0; j < nsides; ++j)
            {
                eid = f_erep[j];
                if (!_old_new_e_map->count(eid))
                    goto NEXT_FACE;
                else
                    eid = _old_new_e_map->find(eid)->second;
                f_erep[j] = eid;
            }

        // only need to handle element *finiteness* when loading entire VD
        if (!_vol)
            if (f_erep[f_erep.size() - 1] == -1)
            {
                // grab adjacent 2 half-infinite edges to create this new edge
                auto e1 = geom().getEdge(f_erep.front());
                auto e2 = geom().getEdge(*(f_erep.end() - 2));
                auto new_e = util::makeEdge(e1[1], e2[1]);
                m_geom.appendEdge(new_e);
                m_is_finite_e.push_back(false);
                f_erep.back() = geom().numEdges() - 1;
                m_is_finite_f.push_back(false);
            }
            else
            {
                m_is_finite_f.push_back(true);
            }

        f_of_vts.resize(nsides);
        trace_face(f_erep, m_geom.m_edges, f_of_vts);

        // add this face to geometry
        // m_geom.appendFace( f_of_vts, f_erep );
        m_geom.appendFace(f_of_vts);
        // save two sites for this face
        m_face_sites.push_back(ivec2(s1, s2));

    NEXT_FACE: // label used to quickly skip processing of cur face
               ;
    }

    m_face_sites.shrink_to_fit();
    m_is_finite_f.shrink_to_fit();
    buffer.clear();
    buffer.shrink_to_fit();
} // VoroInfo::load_voro_faces()

void VoroInfo::infer_sites_from_cell_file(ifstream& _in_cells)
{
    _in_cells.seekg(0, std::ios::end);
    auto len = _in_cells.tellg();
    _in_cells.seekg(std::ios::beg);
    vector<char> buffer(len);
    _in_cells.read(&buffer[0], len);

    int cell_cnt;
    auto delim = " \n\t";
    // read num of cells
    auto token = strtok(&buffer[0], delim);
    cell_cnt = atoi(token);

    vector<int> f_vrep;
    int n_face, f_id;
    vector<point> cell_centroids(cell_cnt);
    for (auto ci = 0; ci < cell_cnt; ++ci)
    {
        token = strtok(nullptr, delim); // cell-id
        token = strtok(nullptr, delim); // num of face index
        n_face = atoi(token);
        int vts_cnt = 0;
        point cent;
        for (auto i = 0; i < n_face; ++i)
        {
            token = strtok(nullptr, delim); // face id
            f_id = atoi(token);
            if (f_id == -1)
                continue;
            geom().getFaceVRep(f_id, f_vrep);
            vts_cnt += f_vrep.size();
            for (auto it = f_vrep.begin(); it != f_vrep.end(); ++it)
            {
                cent += geom().getVert(*it);
            }
        }
        // use cell centroid to find actual sites of voro faces
        if (vts_cnt > 0)
        {
            cent /= vts_cnt;
            cell_centroids[ci] = cent;
        }
        else
        {
            std::cout << "WARNING: couldn't estimate cell centroid for cell "
                      << ci << ". "
                      << "(voro faces' sites will be unstable...)" << std::endl;
        }
    } // end for each cell
    buffer.clear();
    buffer.shrink_to_fit();

    computeInfoRelatedtoSites(cell_centroids);
} // end of load_voro_cells()

void VoroInfo::trace_face(const vector<int>& _f_of_sides,
                          const vector<ivec2>& _edges,
                          vector<int>& _f_of_vts) const
{
    ivec2 e1, e2;
    e1 = _edges[_f_of_sides[0]];
    e2 = _edges[_f_of_sides[1]];
    int pre_v;

    // initialize prev vert
    if (e1[1] != e2[0])
        std::swap(e1[0], e1[1]);
    pre_v = e1[1];
    _f_of_vts[0] = pre_v;

    for (int i = 1; i < _f_of_sides.size(); ++i)
    {
        const auto& e = _edges[_f_of_sides[i]];
        if (e[0] == pre_v)
        {
            pre_v = e[1];
        }
        else
        {
            pre_v = e[0];
        }
        _f_of_vts[i] = pre_v;
    }
}

void VoroInfo::compute_multiplicity_vts()
{
    // construct
    // - edge_id-face_id adjacency table
    vector<vector<int>> e_f_adj_tbl;
    // - vert_id-face_id adjacency table
    vector<vector<int>> v_f_adj_tbl;
    // - count for vert referenced by isolated edges
    vector<int> v_e_refcnt;
    e_f_adj_tbl.resize(m_geom.numEdges());
    v_f_adj_tbl.resize(m_geom.numVts());
    v_e_refcnt.resize(m_geom.numVts(), 0);
    int num_faces = m_geom.numFaces();
    int num_edges = m_geom.numEdges();
    vector<int> f_v_rep, f_e_rep;
    for (int fi = 0; fi < num_faces; ++fi)
    {
        m_geom.getFaceVRep(fi, f_v_rep);
        for (const auto& v : f_v_rep)
        {
            v_f_adj_tbl[v].push_back(fi);
        }
        m_geom.getFaceERep(fi, f_e_rep);
        for (const auto& ei : f_e_rep)
        {
            e_f_adj_tbl[ei].push_back(fi);
        }
    }
    for (size_t ei = 0; ei < num_edges; ++ei)
    {
        if (m_geom.refCount(ei) == 0)
        {
            auto e = m_geom.getEdge(ei);
            v_e_refcnt[e[0]]++;
            v_e_refcnt[e[1]]++;
            // cout << e << " contribute to v_e_refcnt" << endl;
        }
    }

    // track the set of so-far visited elems of cur local star
    set<int> e_visited;
    set<int> f_visited;
    // the queue used to do flooding
    queue<int> flood_queue;

    // the local star ( one-ring faces ) for quick look-up
    set<int> one_ring_faces;

    // identify surface components for each vertex
    m_mult_num_vts.assign(m_geom.numVts(), 1);

    int num_vts = (int)m_geom.numVts();
    for (int vi = 0; vi < num_vts; ++vi)
    {
        if (vi == num_vts - 1)
            [[maybe_unused]] int stop = 1;
        // construct local star for cur vertex
        one_ring_faces.clear();
        one_ring_faces.insert(v_f_adj_tbl[vi].begin(), v_f_adj_tbl[vi].end());

        // start from unvisited face, trace out face components
        e_visited.clear();
        f_visited.clear();
        int ncc = 0;
        for (const auto& fi : one_ring_faces)
        {
            if (f_visited.count(fi) > 0)
                continue;
            ncc++;
            flood_queue.push(fi);
            while (!flood_queue.empty())
            {
                auto cur_fi = flood_queue.front();
                flood_queue.pop();
                if (f_visited.count(cur_fi) > 0)
                    continue;
                f_visited.insert(cur_fi);
                // flood to unvisited neighbor faces, thru unvisited edges
                m_geom.getFaceERep(cur_fi, f_e_rep);
                for (const auto& ei : f_e_rep)
                {
                    if (e_visited.count(ei) == 0)
                    {
                        e_visited.insert(ei);
                        const auto& nb_faces = e_f_adj_tbl[ei];
                        for (const auto& nb_fi : nb_faces)
                        {
                            // only flood to neighbor face (within current
                            // one-ring!)
                            if (f_visited.count(nb_fi) == 0 &&
                                one_ring_faces.count(nb_fi) > 0)
                            {
                                flood_queue.push(nb_fi);
                            }
                        }
                    }
                } // Done: flooding to unvisited neighbor faces
            }
        } // Done: tracing face components

        // now trace out edge components
        // this is simply the number of edges that use cur vertex as one of
        // their ends
        ncc += (ncc == 0 && v_e_refcnt[vi] > 0 ? 1 : 0);

        // now we have traced out all components
        m_mult_num_vts[vi] = ncc;
        // m_mult_num_vts[ vi ] = 1;
        // cout << vi << "'s multiplicity: " << ncc << endl;
        if (m_mult_num_vts[vi] == 0)
            cout << vi << " has no neighborhood" << endl;
    }

    e_f_adj_tbl.clear();
    v_f_adj_tbl.clear();
    e_f_adj_tbl.shrink_to_fit();
    v_f_adj_tbl.shrink_to_fit();
} // compute_multiplicity_vts()

void VoroInfo::recomputeEdgeFlags()
{
    m_edge_valid.resize(m_geom.numEdges());
    m_is_finite_e.resize(m_geom.numEdges());
    for (auto i = 0; i < m_geom.numEdges(); ++i)
    {
        m_edge_valid[i] = computeEdgeValidity(i);
        const auto e = m_geom.getEdge(i);
        m_is_finite_e[i] = isVertexFinite(e[0]) && isVertexFinite(e[1]);
    }
}
void VoroInfo::recomputeFaceFlags()
{
    m_face_valid.resize(m_geom.numFaces());
    m_is_finite_f.resize(m_geom.numFaces());
    vector<int> f_vts;
    for (auto i = 0; i < m_geom.numFaces(); ++i)
    {
        m_face_valid[i] = computeFaceValidity(i);
        m_geom.getFaceVRep(i, f_vts);
        m_is_finite_f[i] = true;
        for (auto j : f_vts)
            if (!isVertexFinite(j))
            {
                m_is_finite_f[i] = false;
                break;
            }
    }
}

} // namespace voxelvoro
