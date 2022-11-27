#include <voxelcore/edgecollapse.h>
#include <voxelcore/geomalgo.h>

#include <algorithm>
#include <cassert>
#include <map>
#include <set>
#include <unordered_set>

namespace voxelvoro
{
using std::map;
using std::set;
using std::unordered_set;
TopoPreservEdgeCollapser::TopoPreservEdgeCollapser(
    const shared_ptr<vector<point>>& _vts_ptr,
    const shared_ptr<vector<uTriFace>>& _tris_ptr)
{
    m_vts_ptr = _vts_ptr;
    m_tris_ptr = _tris_ptr;
    // extract edges
    set<ivec2> edges_set;
    for (const auto& f : *m_tris_ptr)
    {
        edges_set.insert(util::makeEdge(f[0], f[1]));
        edges_set.insert(util::makeEdge(f[1], f[2]));
        edges_set.insert(util::makeEdge(f[0], f[2]));
    }
    m_edges = std::make_shared<vector<ivec2>>();
    m_edges->assign(edges_set.begin(), edges_set.end());
    edges_set.clear();
}
TopoPreservEdgeCollapser::~TopoPreservEdgeCollapser() {}

void TopoPreservEdgeCollapser::preProcess() { compute_order(); }

void TopoPreservEdgeCollapser::collapse()
{
    /*global: construct nbvts, nbfaces for each vertex*/
    // mapping: v id -> list of {nb v id, and nb edge id}
    vector<vector<int>> v_v_adj_tbl;
    v_v_adj_tbl.resize(m_vts_ptr->size());
    // mapping: v id -> list of nb faces ids
    vector<vector<int>> v_f_adj_tbl;
    v_f_adj_tbl.resize(m_vts_ptr->size(), vector<int>());
    // mapping: e id -> list of the faces sharing this edge
    map<ivec2, vector<int>> e_f_adj_tbl;
    typedef map<ivec2, vector<int>>::iterator EFTableIter;
    // face removed tag
    vector<bool> removed_tag_faces;
    removed_tag_faces.resize(m_tris_ptr->size(), false);
    // fill-in those structures
    ivec2 edges_from_face[3];
    for (int fi = 0; fi < m_tris_ptr->size(); ++fi)
    {
        const auto& f = (*m_tris_ptr)[fi];
        util::makeEdgesFromTri(f, edges_from_face);
        for (int i = 0; i < 3; ++i)
        {
            const auto& e = edges_from_face[i];

            if (e == ivec2(727054, 727251))
            {
                std::cout << "edge " << e << " is from face " << f << std::endl;
            }

            // fill-in e-v-adj
            e_f_adj_tbl[e].push_back(fi);
            // fill-in v-f-adj
            v_f_adj_tbl[f[i]].push_back(fi);
        }
    }
    // fill-in v-v-adj
    for (const auto& pr : e_f_adj_tbl)
    {
        const auto& e = pr.first;
        v_v_adj_tbl[e[0]].push_back(e[1]);
        v_v_adj_tbl[e[1]].push_back(e[0]);
    }

    /*collapse all degenerate edges*/
    std::cout << "total # edges to inspect: " << e_f_adj_tbl.size()
              << std::endl;
    auto are_same = [](const point& _v1, const point& _v2) -> bool
    { return util::closeVertices(_v1, _v2, 0.0000001f); };
    bool u_in_bndry[3];
    bool v_in_bndry[3];
    unordered_set<int> lk_u_vts;
    unordered_set<ivec2, ivec2Hash> lk_u_edges;
    ivec2 oppo_e, nb_e;
    unordered_set<int> faces_to_remove;
    vector<int> faces_of_merged_nb;
    unordered_set<int> vts_of_merged_nb;
    int w = -1; // the dummy vertex
    bool edge_collapsed = false;
    // iterate through each edge
    int progress = 0;
    auto update_it_to_next_edge =
        [&e_f_adj_tbl](EFTableIter& it, bool edge_collapsed)
    {
        if (edge_collapsed == true)
        {
            // delete this edge
            it = e_f_adj_tbl.erase(it);
        }
        else
        {
            it++;
        }
    };
    for (auto it = e_f_adj_tbl.begin(); it != e_f_adj_tbl.end(); /*++it*/)
    {
        if (progress++ % 1000000 == 0)
        {
            std::cout << "edges inspected: " << progress << std::endl;
        }
        const auto& e = it->first;
        const auto &u = e[0], &v = e[1];

        // only try to collapse degenerate edge
        if (!are_same((*m_vts_ptr)[u], (*m_vts_ptr)[v]))
        {
            edge_collapsed = false;
            update_it_to_next_edge(it, edge_collapsed);
            continue;
        }

        /* now we know this is a degenerate edge. check link conditions. */
        const auto& lk_uv = it->second;
        const auto& nb_vts_u = v_v_adj_tbl[u];
        const auto& nb_vts_v = v_v_adj_tbl[v];
        const auto& nb_faces_u = v_f_adj_tbl[u];
        const auto& nb_faces_v = v_f_adj_tbl[v];
        assert(!nb_faces_u.empty() && !nb_faces_v.empty());
        for (int bi = 0; bi < 3; ++bi)
        {
            u_in_bndry[bi] = m_vert_ord[u] == bi;
            v_in_bndry[bi] = m_vert_ord[v] == bi;
        }

        /* link condition 1: lk_u intersected with lk_v = lk_uv */
        // check if this holds for the link containing only vertices (nb vts):
        // find the size of the intersected sets, compare the size with that of
        // the edge uv
        lk_u_vts.clear();
        lk_u_vts.insert(nb_vts_u.begin(), nb_vts_u.end());
        int num_common_vts = std::count_if(
            nb_vts_v.begin(), nb_vts_v.end(),
            [&](int _vid) { return lk_u_vts.find(_vid) != lk_u_vts.end(); });
        // if u and v are not manifold vts (i.e. in 1st or higher boundary),
        // then their respective lk includes the dummy vertex w
        num_common_vts += int(!u_in_bndry[0] && !v_in_bndry[0]);
        // if uv is a non-manifold edge, its lk includes the dummy vertex w
        int size_lk_uv = lk_uv.size() + int(lk_uv.size() != 2);
        if (num_common_vts != size_lk_uv)
        {
            edge_collapsed = false;
            update_it_to_next_edge(it, edge_collapsed);
            continue;
        }
        // check if this holds for the link for u, v containing edges (nb
        // edges): reject edge if the intersection between lk_u and lk_v is
        // non-empty (as lk_uv contains only vertices on 2-complex)
        lk_u_edges.clear();
        bool common_edge_non_empty = false;
        for (const auto& nb_fi : nb_faces_u)
        {
            oppo_e = util::oppoEdge((*m_tris_ptr)[nb_fi], u);
            lk_u_edges.insert(oppo_e);
        }
        for (const auto& nb_fi : nb_faces_v)
        {
            oppo_e = util::oppoEdge((*m_tris_ptr)[nb_fi], v);
            if (lk_u_edges.count(oppo_e) > 0)
            {
                // found common edge. so current edge is not collapsible
                common_edge_non_empty = true;
                break;
            }
        }
        if (common_edge_non_empty == false)
        {
            // also account for the common edges with dummy vertex as one of the
            // ends for u, v find the edges involving dummy vertex in their
            // resp. lk and see if there are any in common find such edges for u
            // first
            for (const auto& nb_vi : nb_vts_u)
            {
                nb_e = util::makeEdge(u, nb_vi);
                // if nb edge is singular, a dummy vertex is connected to the
                // nb_vi to form an oppo_e for vertex u. mark this oppo-e as a
                // potential common edge between lk_u and lk_v
                if (e_f_adj_tbl.find(nb_e)->second.size() != 2)
                {
                    lk_u_edges.insert(util::makeEdge(nb_vi, w));
                }
            }
            // do similar thing to vertex v:
            // find its dummy-connected lk edges and check each one of them
            // to see if there are any common edges with lk_u
            for (const auto& nb_vi : nb_vts_v)
            {
                nb_e = util::makeEdge(v, nb_vi);
                // if nb edge is singular, a dummy vertex is connected to the
                // nb_vi to form an oppo_e for vertex v. check if lk_u has any
                // edge in common
                if (e_f_adj_tbl.find(nb_e)->second.size() != 2)
                {
                    if (lk_u_edges.count(util::makeEdge(nb_vi, w)) > 0)
                    {
                        common_edge_non_empty = true;
                        break;
                    }
                }
            }
        }
        if (common_edge_non_empty)
        {
            edge_collapsed = false;
            update_it_to_next_edge(it, edge_collapsed);
            continue;
        }

        // link condition 2:
        // lk_1_u intersected with lk_1_v == empty set
        bool lk_1_common_non_empty = false;
        lk_u_vts.clear();
        for (const auto& nb_vi : nb_vts_u)
        {
            nb_e = util::makeEdge(u, nb_vi);
            // if nb edge is singular, put the nb to lk_1_u
            if (e_f_adj_tbl.find(nb_e)->second.size() != 2)
            {
                lk_u_vts.insert(nb_vi);
            }
        }
        for (const auto& nb_vi : nb_vts_v)
        {
            nb_e = util::makeEdge(v, nb_vi);
            // if nb edge is singular, check the nb to see if it's common in
            // lk_1_u
            if (e_f_adj_tbl.find(nb_e)->second.size() != 2)
            {
                if (lk_u_vts.count(nb_vi) > 0)
                {
                    lk_1_common_non_empty = true;
                    break;
                }
            }
        }
        // also need to check if dummy w is common to lk1_u, and lk1_v
        if (u_in_bndry[2] && v_in_bndry[2])
            lk_1_common_non_empty = true;
        // if the common part is not empty, reject the edge
        if (lk_1_common_non_empty)
        {
            edge_collapsed = false;
            update_it_to_next_edge(it, edge_collapsed);
            continue;
        }

        /*the edge uv has passed all tests, meaning we should collapse it:
        collapse the edge. merge/update the affected vts/edges' neighborhoods
        properly including the v-v-adj, v-f-adj, and e-v-adj */
        // for each vert in lk_uv, remove its nb faces that are "to-remove"
        // put those faces in a set for quick query later
        faces_to_remove.clear();
        for (const auto& fi_to_remove : lk_uv)
        {
            const auto& f = (*m_tris_ptr)[fi_to_remove];
            const auto& s = util::oppoVertex(f, e);
            faces_to_remove.insert(fi_to_remove);
            removed_tag_faces[fi_to_remove] = true;
            // remove "removed" faces from nb faces of vert in lk_uv
            auto& nb_faces_s = v_f_adj_tbl[s];
            for (auto it = nb_faces_s.begin(); it != nb_faces_s.end(); ++it)
            {
                if (*it == fi_to_remove)
                {
                    nb_faces_s.erase(it);
                    break;
                }
            }
            // merge/update nb faces for edges us, and vs
            // assign merged faces to the nbfaces of edge vs
            faces_of_merged_nb.clear();
            ivec2 e_us = util::makeEdge(u, s);
            auto e_us_it = e_f_adj_tbl.find(e_us);
            if (e_us_it == e_f_adj_tbl.end())
            {
                std::cout << "Error (during edge u-s lookup): edge "
                          << e_us_it->first << " is not in e-f-adj-tbl."
                          << std::endl;
                exit(-1);
            }
            const auto& nb_faces_us = e_us_it->second;
            for (auto it = nb_faces_us.begin(); it != nb_faces_us.end(); ++it)
            {
                if (*it != fi_to_remove)
                {
                    faces_of_merged_nb.push_back(*it);
                }
            }
            ivec2 e_vs = util::makeEdge(v, s);
            auto e_vs_it = e_f_adj_tbl.find(e_vs);
            if (e_vs_it == e_f_adj_tbl.end())
            {
                std::cout << "Error (during edge v-s lookup): edge "
                          << e_vs_it->first << " is not in e-f-adj-tbl."
                          << std::endl;
                exit(-1);
            }
            const auto& nb_faces_vs = e_vs_it->second;
            for (auto it = nb_faces_vs.begin(); it != nb_faces_vs.end(); ++it)
            {
                if (*it != fi_to_remove)
                {
                    faces_of_merged_nb.push_back(*it);
                }
            }
            e_vs_it->second = faces_of_merged_nb;

            // remove "u" from s' nbvts
            auto& nbvts_s_to_modify = v_v_adj_tbl[s];
            for (auto it = nbvts_s_to_modify.begin();
                 it != nbvts_s_to_modify.end(); ++it)
            {
                if (*it == u)
                {
                    nbvts_s_to_modify.erase(it);
                    break;
                }
            }
        }

        // rename all occurrences of u to v in the nb faces of u
        for (const auto fi : nb_faces_u)
        {
            auto& f = (*m_tris_ptr)[fi];
            util::replaceVinF(f, u, v);
        }
        // merge u's nbvts (excluding v) to v's nbvts
        vts_of_merged_nb.clear();
        auto& nbvts_v_to_modify = v_v_adj_tbl[v];
        for (const auto& nbvi : nb_vts_u)
        {
            if (nbvi != v)
                vts_of_merged_nb.insert(nbvi);
        }
        for (const auto& nbvi : nb_vts_v)
        {
            if (nbvi != u)
                vts_of_merged_nb.insert(nbvi);
        }
        nbvts_v_to_modify.clear();
        nbvts_v_to_modify.insert(nbvts_v_to_modify.end(),
                                 vts_of_merged_nb.begin(),
                                 vts_of_merged_nb.end());
        // also need to replace u with v for each nb vert of u
        // and replace old nb edge nb-u with nb-v in the e-f-adj-tbl,
        // where nb is u's nb (excluding v)
        for (const auto& nbvi : nb_vts_u)
        {
            auto& nbs = v_v_adj_tbl[nbvi];
            for (auto& vi : nbs)
            {
                if (vi == u)
                {
                    vi = v;
                    break;
                }
            }
            // replace
            if (nbvi != v)
            {
                auto old_e_it = e_f_adj_tbl.find(util::makeEdge(nbvi, u));
                ivec2 new_e = util::makeEdge(nbvi, v);
                auto new_e_it = e_f_adj_tbl.find(new_e);
                if (new_e_it == e_f_adj_tbl.end())
                    e_f_adj_tbl[new_e] = old_e_it->second;
                if (old_e_it == e_f_adj_tbl.end())
                {
                    std::cout << "Error (during edge u_nb-u lookup): edge "
                              << old_e_it->first << " is not in e-f-adj-tbl."
                              << std::endl;
                    exit(-1);
                }
                e_f_adj_tbl.erase(old_e_it);
            }
        }
        v_v_adj_tbl[u].clear();

        // merge the nb faces of u into v (excluding faces-to-remove)
        faces_of_merged_nb.clear();
        for (const auto fi : nb_faces_u)
        {
            if (faces_to_remove.count(fi) == 0)
            {
                faces_of_merged_nb.push_back(fi);
            }
        }
        for (const auto fi : nb_faces_v)
        {
            if (faces_to_remove.count(fi) == 0)
            {
                faces_of_merged_nb.push_back(fi);
            }
        }
        v_f_adj_tbl[u].clear();
        v_f_adj_tbl[v] = faces_of_merged_nb;

        // after merging, update v's order
        m_vert_ord[v] = std::max(m_vert_ord[u], m_vert_ord[v]);

        edge_collapsed = true;

        update_it_to_next_edge(it, edge_collapsed);
    } // for () : try to collapse all degenerate edges

    std::cout << "removing faces ..." << std::endl;

    /*now actually remove the faces marked as "removed"*/
    vector<uTriFace> remaining_faces;
    for (int fi = 0; fi < m_tris_ptr->size(); ++fi)
    {
        if (removed_tag_faces[fi] == false)
        {
            remaining_faces.push_back((*m_tris_ptr)[fi]);
        }
    }
    *m_tris_ptr = remaining_faces;
} // TopoPreservEdgeCollapser::collapse()

// compute topo order by constructing simple version of topo graph locally
void TopoPreservEdgeCollapser::compute_order()
{
    /*one-ring structures needed for topo graph computation*/
    // mapping: global v id -> local one-ring v id
    map<int, int> v_map;
    // local one-ring vts id -> map global vts id
    vector<int> vts;
    vts.reserve(100);
    // one-ring edges
    vector<ivec2> edges;
    edges.reserve(600);
    // edge visited tag
    vector<bool> visited_e;
    visited_e.reserve(edges.capacity());
    // one-ring faces
    vector<uTriFace> faces;
    faces.reserve(200);
    // one-ring map: v -> list of {nb u and edge id for <v, u>}
    vector<vector<ivec2>> v_v_adj_tbl;
    v_v_adj_tbl.reserve(vts.capacity());
    // remaining degree count for each vert in one-ring
    // - degree reduces by one when a nb edge is visited;
    // - degree reduces to zero when a vertex is fully visited;
    vector<int> v_remain_deg;
    v_remain_deg.reserve(vts.capacity());

    /*one-ring topo structures*/
    // topo arcs for the one-ring neighborhood
    vector<ivec2> topo_edges;
    topo_edges.reserve(100);
    // topo vts that topo arcs share
    vector<int> topo_vts;
    topo_vts.reserve(100);
    // map: topo v -> list of {nb topo u and topo edge id for <v, u>}
    vector<vector<ivec2>> topo_v_v_adj_tbl;
    topo_v_v_adj_tbl.reserve(vts.capacity());
    // visited tag for topo vert
    vector<bool> visited_topo_v;
    visited_topo_v.reserve(vts.capacity());

    /* global structures */
    // map: v -> list of nb face ids
    vector<vector<int>> all_v_face_adj_tbl;
    // degree count for edge
    map<ivec2, int> all_e_degree_cnt;

    /*initialization work*/
    // collect nb faces for each vert
    // and ref cnt for each edge
    all_v_face_adj_tbl.resize(m_vts_ptr->size());
    for (auto fi = 0; fi < m_tris_ptr->size(); ++fi)
    {
        const auto& f = (*m_tris_ptr)[fi];
        all_v_face_adj_tbl[f[0]].push_back(fi);
        all_v_face_adj_tbl[f[1]].push_back(fi);
        all_v_face_adj_tbl[f[2]].push_back(fi);
        all_e_degree_cnt[util::makeEdge(f[0], f[1])]++;
        all_e_degree_cnt[util::makeEdge(f[1], f[2])]++;
        all_e_degree_cnt[util::makeEdge(f[0], f[2])]++;
    }

    /*iteration stage*/
    // determine the order for each vertex
    for (auto vi = 0; vi < m_vts_ptr->size(); ++vi)
    {
        /* construct one-ring quantities */
        v_map.clear();
        vts.clear();
        edges.clear();
        faces.clear();
        v_remain_deg.clear();
        // pick out of global vts list the one-ring vts for current vertex vi
        // also construct one-ring edges between them properly
        size_t pre_size = 0;
        const auto& nb_faces = all_v_face_adj_tbl[vi];
        ivec2 oppo_e;
        int oppo_v1, oppo_v2;
        for (const auto& fi : nb_faces)
        {
            const auto& f = (*m_tris_ptr)[fi];
            for (int i = 0; i < 3; ++i)
            {
                if (f[i] == vi)
                {
                    oppo_v1 = f[(i + 1) % 3];
                    oppo_v2 = f[(i + 2) % 3];
                    pre_size = v_map.size();
                    auto& v1_local_id = v_map[oppo_v1];
                    if (v_map.size() > pre_size)
                    {
                        vts.push_back(oppo_v1);
                        v1_local_id = vts.size() - 1;
                    }
                    pre_size = v_map.size();
                    auto& v2_local_id = v_map[oppo_v2];
                    if (v_map.size() > pre_size)
                    {
                        vts.push_back(oppo_v2);
                        v2_local_id = vts.size() - 1;
                    }

                    // also, add this oppo edge to the v-v adj table
                    oppo_e = util::makeEdge(v1_local_id, v2_local_id);
                    edges.push_back(oppo_e);
                    break;
                }
            }
        }

        // build adj table for one-ring vts
        v_v_adj_tbl.clear();
        v_v_adj_tbl.resize(vts.size());
        visited_e.clear();
        for (int i = 0; i < edges.size(); ++i)
        {
            const auto& e = edges[i];
            visited_e.push_back(false);
            v_v_adj_tbl[e[0]].push_back(ivec2(e[1], i));
            v_v_adj_tbl[e[1]].push_back(ivec2(e[0], i));
        }

        auto is_nonmanifold_v = [&](int _i) -> bool
        { return v_v_adj_tbl[_i].size() != 2; };
        bool has_nonMani_vert = false;

        // set degree cnt for one-ring vts
        for (int i = 0; i < vts.size(); ++i)
        {
            v_remain_deg.push_back(v_v_adj_tbl[i].size());
            if (v_v_adj_tbl[i].size() != 2)
            {
                has_nonMani_vert = true;
            }
        }

        /* now build topo graph:
        - trace topo arcs from non-manifold vert to non-manifold vert
        - build topo vts adj table
        */
        topo_vts.clear();
        topo_edges.clear();
        topo_v_v_adj_tbl.clear();
        topo_v_v_adj_tbl.resize(vts.size(), vector<ivec2>());
        int start_v, cur_v, next_v, end_v;
        for (auto i = 0; i < vts.size(); ++i)
        {
            // set topo vts as junction(non-manifold) vts
            if (is_nonmanifold_v(i))
                topo_vts.push_back(i);
        }
        for (auto i = 0; i < vts.size(); ++i)
        {
            // if no non-manifold vert in one-ring,
            // then any unvisited vert can serve as start point
            // otherwise, we only start from an unvisited non-manifold vert
            if (!has_nonMani_vert)
            {
                start_v = i;
            }
            else if (is_nonmanifold_v(i))
            {
                start_v = i;
            }
            else
            {
                continue;
            }

            while (v_remain_deg[start_v] > 0)
            {
                // initialize tracing: start from start_v to trace
                next_v = cur_v = start_v;
                // by following the unvisited nb edge thru manifold vert,
                // tracing stops when either next vert is non-manifold,
                // or next vert's remaining deg is 0
                while (v_remain_deg[next_v] > 0)
                {
                    cur_v = next_v;
                    const auto& nbs = v_v_adj_tbl[cur_v];
                    for (const auto& nb : nbs)
                    {
                        if (!visited_e[nb[1]])
                        {
                            next_v = nb[0];
                            visited_e[nb[1]] = true;
                            v_remain_deg[cur_v]--;
                            v_remain_deg[next_v]--;
                            break;
                        }
                    }

                    if (is_nonmanifold_v(next_v))
                        break;
                }
                end_v = next_v;

                // now we get the start and end points of a topo arc
                ivec2 topo_e;
                if (start_v == end_v) // a loop
                {
                    if (is_nonmanifold_v(start_v)) // regular loop
                    {
                        topo_e = ivec2(start_v, start_v);
                    }
                    else // pure loop (loop w/o ends)
                    {
                        topo_e = ivec2(-1, -1);
                    }
                }
                else // regular topo edge
                {
                    topo_e = util::makeEdge(start_v, end_v);
                }
                topo_edges.push_back(topo_e);

                // add non-loop adj info to topo adj tbl
                if (topo_e[0] != topo_e[1])
                {
                    topo_v_v_adj_tbl[topo_e[0]].push_back(
                        ivec2(topo_e[1], topo_edges.size() - 1));
                    topo_v_v_adj_tbl[topo_e[1]].push_back(
                        ivec2(topo_e[0], topo_edges.size() - 1));
                }
            }
        }

        /*now topo graph for this vertex is ready. identify its type/order*/
        // a manifold vert has order 0
        if (topo_edges.size() == 1 && topo_edges[0][0] == -1)
            m_vert_ord.push_back(0);
        else
        {
            bool has_cycle =
                false; // including loop, but can be higher order cycle
            bool has_loop = false; // pure or regular loop of one topo edge
            // check for loop edge first
            for (const auto& te : topo_edges)
            {
                if (te[0] == te[1])
                {
                    has_loop = true;
                    has_cycle = true;
                    break;
                }
            }
            if (!has_cycle)
            {
                // check for higher order cycle
                // Since we know at this point there is no loop topo edge in the
                // graph and all topo vts are non-manifold vts so we can detect
                // cycle easily locally at each topo edge: we set the ends of a
                // topo edge as visited after inspecting that edge if at some
                // point we find both ends of an edge are "visited", then there
                // is a cycle.
                visited_topo_v.assign(vts.size(), false);
                for (const auto& te : topo_edges)
                {
                    if (visited_topo_v[te[0]] && visited_topo_v[te[1]])
                    {
                        has_cycle = true;
                        break;
                    }
                    visited_topo_v[te[0]] = true;
                    visited_topo_v[te[1]] = true;
                }
            }
            // a boundary vert has order 1
            if (!has_cycle)
            {
                m_vert_ord.push_back(1);
            }
            else
            {
                // a vert on non-manifold curve has order 1
                // at this point, there is no cycle,
                // so it suffices to check if there are two topo vts
                // to determine whether it's on non-manifold curve
                if (topo_vts.size() == 2 && !has_loop)
                {
                    m_vert_ord.push_back(1);
                }
                else // it's a junction point that has most complicated
                     // neighborhood
                {
                    m_vert_ord.push_back(2);
                }
            }
        }
    }

    // compute order for each edge
    m_edge_ord.resize(m_edges->size());
    for (int ei = 0; ei < m_edges->size(); ++ei)
    {
        const auto& e = (*m_edges)[ei];
        if (all_e_degree_cnt[e] != 2)
        {
            m_edge_ord[ei] = 1;
        }
        else
        {
            m_edge_ord[ei] = 0;
        }
    }
} // TopoPreservEdgeCollapser::compute_order()

} // namespace voxelvoro
