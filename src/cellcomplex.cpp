#include <queue>
#include <set>
#include <unordered_map>
using std::pair;
using std::queue;
using std::set;
using std::unordered_map;
#include "cellcomplex.h"
#include <TriMesh.h>

cellcomplex::cellcomplex() { invalidateStates(); }

cellcomplex::cellcomplex(const vector<point>& _vts, const vector<ivec2>& _edges,
                         const vector<uTriFace>& _faces,
                         bool _perform_finalize /*= true*/)
{
    invalidateStates();
    for (const auto& v : _vts)
        appendVert(v);
    for (const auto& e : _edges)
        appendEdge(e);
    for (const auto& t : _faces)
        appendFace(t);
    if (_perform_finalize)
        finalize();
}

cellcomplex::~cellcomplex() {}

void cellcomplex::invalidateStates()
{
    m_has_infinite_v = m_has_infinite_e = false;
    m_is_finalized = false;
}

bool cellcomplex::needVVAdjacency()
{
    if (m_nb_edges_for_v.size() == numVts())
        return false;
    // build v-v-adjacency
    init_V_E_adjacency();
    for (size_t ei = 0; ei < numEdges(); ++ei)
    {
        const auto& e = getEdge(ei);
        add_nb_edge(e[0], ei);
        add_nb_edge(e[1], ei);
    }
    return true;
}
bool cellcomplex::needEdges()
{
    if (isFinalized())
        return false;

    vector<int> f_vts;
    std::unordered_set<ivec2, ivec2Hash> edges_set;
    for (auto fi = 0; fi < numFaces(); ++fi)
    {
        f_vts.clear();
        getFaceVRep(fi, f_vts);
        for (auto j = 0; j < f_vts.size(); ++j)
        {
            auto e = util::makeEdge(f_vts[j], f_vts[(j + 1) % f_vts.size()]);
            edges_set.insert(e);
        }
    }
    for (auto e : edges_set)
        appendEdge(e);
    edges_set.clear();

    return true;
}

int cellcomplex::compute_conn_cmpnts(const vector<int>* _subset_faces) const
{
    int ncc = 0;
    auto& ccg = *this;
    // construct
    // - edge_id-face_id adjacency table
    vector<vector<int>> e_f_adj_tbl;
    e_f_adj_tbl.resize(ccg.numEdges());
    vector<bool> visited_f(ccg.numFaces(), false);
    vector<bool> visited_e(ccg.numEdges(), false);
    vector<int> f_e_rep;
    if (_subset_faces)
    {
        // conn cmpnts computation is constrained within only a subset of faces
        // therefore we need to *block* other faces properly
        set<int> subset(_subset_faces->begin(), _subset_faces->end());
        for (int fi = 0; fi < ccg.numFaces(); ++fi)
        {
            if (subset.count(fi) > 0)
            {
                visited_f[fi] = false; // this face is visitable!
                ccg.getFaceERep(fi, f_e_rep);
                for (const auto& ei : f_e_rep)
                {
                    e_f_adj_tbl[ei].push_back(fi);
                }
            }
            else
            {
                visited_f[fi] = true; // block this face, i.e. not visitable!
            }
        }
        subset.clear();
    }
    else
    {
        // computation occurs in the set of all faces
        for (int fi = 0; fi < ccg.numFaces(); ++fi)
        {
            visited_f[fi] = false;
            ccg.getFaceERep(fi, f_e_rep);
            for (const auto& ei : f_e_rep)
            {
                e_f_adj_tbl[ei].push_back(fi);
            }
        }
    }

    // trace components starting from unvisited face
    queue<int> face_q;
    for (int fi = 0; fi < ccg.numFaces(); ++fi)
    {
        if (visited_f[fi])
            continue;

        ncc++;
        face_q.push(fi);
        while (!face_q.empty())
        {
            auto cur_fi = face_q.front();
            face_q.pop();
            if (visited_f[cur_fi])
                continue;

            visited_f[cur_fi] = true;
            ccg.getFaceERep(cur_fi, f_e_rep);
            for (const auto& ei : f_e_rep)
            {
                if (visited_e[ei])
                    continue;

                visited_e[ei] = true;
                const auto& nb_faces = e_f_adj_tbl[ei];
                for (const auto& nb_fi : nb_faces)
                {
                    if (visited_f[nb_fi])
                        continue;

                    face_q.push(nb_fi);
                }
            }
        }
    }
    visited_f.clear();
    visited_e.clear();
    visited_f.shrink_to_fit();
    visited_e.shrink_to_fit();
    e_f_adj_tbl.clear();
    e_f_adj_tbl.shrink_to_fit();

    return ncc;
}
int cellcomplex::compute_conn_cmpnts(const vector<int>* _subset_faces,
                                     int overload_func_tag) const
{
    int ncc = 0;
    auto& ccg = *this;
    // construct
    // - vert_id-vert_id adjacency table
    vector<vector<int>> v_v_adj_tbl;
    v_v_adj_tbl.resize(ccg.numVts());
    vector<bool> visited(ccg.numVts(), false);
    vector<int> f_v_rep;
    vector<int> f_e_rep;
    if (_subset_faces)
    {
        /*flooding will occur in a subset of vts*/
        // grab all relevant vts & build v-v-adj
        set<int> vts_set;
        for (const auto& fi : *_subset_faces)
        {
            ccg.getFaceVRep(fi, f_v_rep);
            vts_set.insert(f_v_rep.begin(), f_v_rep.end());
            // add v-v adjacency
            for (auto j = 0; j < f_v_rep.size(); ++j)
            {
                auto e = util::makeEdge(f_v_rep[j],
                                        f_v_rep[(j + 1) % f_v_rep.size()]);
                v_v_adj_tbl[e[0]].push_back(e[1]);
                v_v_adj_tbl[e[1]].push_back(e[0]);
            }
        }
        // properly set other vts to has-been-visited before flooding
        // to make them *non-visitable*,
        // since our flooding is constrained on a subset of vts
        for (int vi = 0; vi < ccg.numVts(); ++vi)
        {
            if (vts_set.count(vi) == 0)
                visited[vi] = true;
        }
        vts_set.clear();
    }
    else
    {
        /* flooding will consider all vts */
        // build v-v-adj table
        for (const auto& e : ccg.m_edges)
        {
            v_v_adj_tbl[e[0]].push_back(e[1]);
            v_v_adj_tbl[e[1]].push_back(e[0]);
        }
    }

    ncc = util::findNumConnComponents(v_v_adj_tbl, ccg.numVts(), visited);
    return ncc;
}

void cellcomplex::appendFace(const vector<int>& _f_of_vts)
{
    int nsides = _f_of_vts.size();
    m_faceVtsN.push_back(nsides);
    auto& f_vrep = m_faces_v_rep[nsides];
    // save the index of the start of this face
    m_faceIndex.push_back(f_vrep.size());
    // add to v-rep
    f_vrep.insert(f_vrep.end(), _f_of_vts.begin(), _f_of_vts.end());
}

void cellcomplex::appendFace(const uTriFace& _f_tri)
{
    int nsides = 3;
    m_faceVtsN.push_back(nsides);
    auto& f_vrep = m_faces_v_rep[nsides];
    m_faceIndex.push_back(f_vrep.size());
    f_vrep.insert(f_vrep.end(), _f_tri.begin(), _f_tri.end());
}

void cellcomplex::appendFace(const vector<int>& _f_of_vts,
                             const vector<int>& _f_of_sides)
{
    int nsides = _f_of_vts.size();
    m_faceVtsN.push_back(nsides);
    auto& f_vrep = m_faces_v_rep[nsides];
    auto& f_erep = m_faces_e_rep[nsides];
    // save the index of the start of this face
    m_faceIndex.push_back(f_vrep.size());
    // add to v-rep & e-rep
    f_vrep.insert(f_vrep.end(), _f_of_vts.begin(), _f_of_vts.end());
    f_erep.insert(f_erep.end(), _f_of_sides.begin(), _f_of_sides.end());
}

size_t cellcomplex::createInfiniteVertex()
{
    if (m_has_infinite_v)
        return m_inf_v_id;
    auto bbox = trimesh::box();
    for (const auto& p : m_vts)
        bbox += p;
    auto inf_vertex = bbox.center() + (bbox.max - bbox.min);
    m_vts.push_back(inf_vertex);
    m_has_infinite_v = true;
    m_inf_v_id = m_vts.size() - 1;
    return m_inf_v_id;
}

size_t cellcomplex::createInfiniteEdge()
{
    if (m_has_infinite_e)
        return m_inf_e_id;
    auto inf_v_id = createInfiniteVertex();
    m_edges.push_back(util::makeEdge(inf_v_id, inf_v_id));
    m_has_infinite_e = true;
    m_inf_e_id = m_edges.size() - 1;
    return m_inf_e_id;
}

size_t cellcomplex::createInfiniteVertex(const point& _p)
{
    m_vts.push_back(_p);
    m_inf_v_id = m_vts.size() - 1;
    m_has_infinite_v = true;
    return m_inf_v_id;
}

void cellcomplex::clear()
{
    invalidateStates();
    clearVertList();
    clearEdgeList();
    clearFaceList();
    clear_adjacency();
}
void cellcomplex::clearFaceList()
{
    invalidateStates();
    m_faceIndex.clear();
    m_faceVtsN.clear();
    m_faces_v_rep.clear();
    m_faces_e_rep.clear();
}
void cellcomplex::clearEdgeList()
{
    invalidateStates();
    m_edges.clear();
}
void cellcomplex::clearVertList()
{
    invalidateStates();
    m_vts.clear();
}

vector<int> cellcomplex::refCntPerVert() const
{
    vector<int> ref_cnt(numVts(), 0);
    for (auto i = 0; i < ref_cnt.size(); ++i)
    {
        ref_cnt[i] = cntNbEdgesofVert(i);
    }
    return ref_cnt;
}
vector<int> cellcomplex::refCntPerEdge() const
{
    vector<int> ref_cnt(numEdges(), 0);
    for (auto i = 0; i < ref_cnt.size(); ++i)
    {
        ref_cnt[i] = cntNbFacesofEdge(i);
    }
    return ref_cnt;
}

void cellcomplex::init_V_E_adjacency()
{
    m_nb_edges_for_v.assign(m_vts.size(), vector<int>());
}
void cellcomplex::init_E_F_adjacency()
{
    m_nb_faces_for_e.assign(m_edges.size(), vector<int>());
}
void cellcomplex::init_adjacency()
{
    init_V_E_adjacency();
    init_E_F_adjacency();
}

void cellcomplex::add_nb_edge(int _vi, int _ei)
{
    m_nb_edges_for_v[_vi].push_back(_ei);
}

void cellcomplex::add_nb_face(int _ei, int _nb_fi)
{
    m_nb_faces_for_e[_ei].push_back(_nb_fi);
}

void cellcomplex::clear_adjacency()
{
    m_nb_edges_for_v.clear();
    m_nb_faces_for_e.clear();
}

void cellcomplex::finalize()
{
    if (isFinalized())
        return;
    /* time stats */
    timer t_total;
    t_total.start();

    timer init_edges_t, extract_edges_t, save_edges_t;
    timer init_adj_t, build_adj_t;
    // check if the edges info has been finalized or not
    if (m_faces_e_rep.empty())
    { // need to finalize
        /*auto e_debug = geom::makeEdge( 1, 60 );
        vector<int> f_vts_debug = { 33,1,60,35 };*/
        // complete edge-idx rep. for faces
        unordered_map<ivec2, int, ivec2Hash> edge_index_map;
        vector<int> vts_of_f;
        vector<ivec2> edges_of_f;

        // init edge-index-map with existing edges before we add new edges
        // notice: order of edges will be kept as index
        init_edges_t.start();
        auto prev_num_of_edges = m_edges.size();
        for (size_t ei = 0; ei < m_edges.size(); ++ei)
        {
            auto pre_size = edge_index_map.size();
            auto& e_idx = edge_index_map[m_edges[ei]];
            auto after_size = edge_index_map.size();
            if (after_size > pre_size)
            {
                e_idx = pre_size;
            }
        }
        init_edges_t.stop();

        // extract edges used by faces
        extract_edges_t.start();
        for (size_t fi = 0; fi < numFaces(); ++fi)
        {
            // we need to extract all edges of faces
            getFaceVRep(fi, vts_of_f);
            /*if ( vts_of_f == f_vts_debug )
                    int stop = 1;*/
            util::traceFace(vts_of_f, edges_of_f);
            // save the edge-id rep. of the face into face-list
            for (size_t ei = 0; ei < edges_of_f.size(); ++ei)
            {
                const auto& e = edges_of_f[ei];
                /*if ( e == e_debug )
                        int stop = 1;*/
                int e_idx;
                auto find_it = edge_index_map.find(e);
                if (find_it != edge_index_map.end())
                    e_idx = find_it->second;
                else
                { // this is a new edge. remember its index.
                    e_idx = edge_index_map.size();
                    edge_index_map[e] = e_idx;
                }
                // actually save the edge-id
                m_faces_e_rep[edges_of_f.size()].push_back(e_idx);
            }
        }
        extract_edges_t.stop();

        // report: are there new edges extracted from face sides, in addition to
        // the original edges?
        if (m_edges.size() < edge_index_map.size())
            std::cout << "INFO: # face sides about to add as new edges = "
                      << edge_index_map.size() - m_edges.size() << endl;

        // save each edge at its corresponding index in edge-list
        save_edges_t.start();
        m_edges.resize(edge_index_map.size());
        for (const auto& e_id_pair : edge_index_map)
        {
            /*if ( e_id_pair.first == e_debug )
                    int stop = 1;*/
            m_edges[e_id_pair.second] = e_id_pair.first;
        }
        edge_index_map.clear();
        save_edges_t.stop();
    }

    // ready to build adjacency info
    init_adj_t.start();
    init_adjacency();
    init_adj_t.stop();

    // update adjacency info for vts
    build_adj_t.start();
    for (size_t ei = 0; ei < numEdges(); ++ei)
    {
        const auto& e = getEdge(ei);
        add_nb_edge(e[0], ei);
        add_nb_edge(e[1], ei);
    }

    // update adjacency info for edges
    m_edge_ref_cnt.resize(m_edges.size(), 0);
    vector<int> e_indices_of_f;
    for (size_t fi = 0; fi < numFaces(); ++fi)
    {
        getFaceERep(fi, e_indices_of_f);
        for (const auto& eid : e_indices_of_f)
        {
            m_edge_ref_cnt[eid]++; // # ref by faces
            add_nb_face(eid, fi);  // nb face add to the edge's list
        }
    }
    build_adj_t.stop();
    t_total.stop();

    /*cout << "init edge struct: "
            << init_edges_t.elapseMilli().count() << "ms" << endl;
    cout << "extract edges from faces: "
            << extract_edges_t.elapseMilli().count() << "ms" << endl;
    cout << "store edges in order: "
            << save_edges_t.elapseMilli().count() << "ms" << endl;
    cout << "init adj list: "
            << init_adj_t.elapseMilli().count() << "ms" << endl;
    cout << "build adj list: "
            << build_adj_t.elapseMilli().count() << "ms" << endl;*/
    cout << "time -> finalize CC: " << t_total.elapseMilli().count() << " ms"
         << endl;
    m_is_finalized = true;
}

void cellcomplex::eulerChar(eulerchar& _ec) const
{
    int ncc = util::findNumConnComponents(this->m_vts, this->m_edges);
    _ec.V = numVts();
    _ec.E = numEdges();
    _ec.F = numFaces();
    _ec.T = 0;
    _ec.euler = _ec.V - _ec.E + _ec.F - _ec.T;
    _ec.C = ncc;
}
void cellcomplex::eulerChar(eulerchar& _ec)
{
    needEdges();
    int ncc = util::findNumConnComponents(this->m_vts, this->m_edges);
    _ec.V = numVts();
    _ec.E = numEdges();
    _ec.F = numFaces();
    _ec.T = 0;
    _ec.euler = _ec.V - _ec.E + _ec.F - _ec.T;
    _ec.C = ncc;
}