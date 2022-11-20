#include <filesystem>
//#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <trimesh/TriMesh.h>

#include <voxelcore/exporters.h>
#include <voxelcore/graphapp.h>
#include <voxelcore/plyall.h>

namespace fs = std::filesystem;

namespace graphapp
{
bool graphapp::readGraph(const string& _filename, vector<point>& _nodes,
                         vector<ivec2>& _edges, vector<vector<float>>& _msures)
{
    // verify if format supported or not
    if (fs::path(_filename).extension() != ".skMsure")
        return false;

    std::ifstream ifile(_filename);
    if (!ifile)
        return false;
    int dump;
    int nNodes, nEdges, nFaces; // faces are not stored now
    ifile >> nNodes >> nEdges >> nFaces;

    _nodes.resize(nNodes);
    point p;
    for (auto i = 0; i < nNodes; ++i)
    {
        ifile >> p[0] >> p[1] >> p[2];
        _nodes[i] = p;
    }

    _edges.clear();
    _edges.reserve(nEdges);
    _msures.resize(3);
    auto& bt3 = _msures[0];
    auto& bt2 = _msures[1];
    auto& bt1 = _msures[2];
    bt3.resize(nEdges);
    bt2.resize(nEdges);
    bt1.resize(nEdges);
    ivec2 e;
    for (auto i = 0; i < nEdges; ++i)
    {
        ifile >> e[0] >> e[1];
        _edges.push_back(e);
        ifile >> bt3[i] >> bt2[i] >> bt1[i];
    }

    // exhaust faces
    for (auto i = 0; i < nFaces; ++i)
    {
        // 3 indices and 2 measure values
        ifile >> dump >> dump >> dump >> dump >> dump;
    }

    return true;
}

bool readGraph(const string& _filename, WeightedGraph& _bg,
               vector<vector<float>>& _e_msures)
{
    vector<point> nodes;
    vector<ivec2> edges;
    auto suc = readGraph(_filename, nodes, edges, _e_msures);
    if (!suc)
        return false;

    vector<float> wts;
    computeWeights(nodes, edges, _e_msures, Weight::BT3, wts);
    makeBoostGraph(nodes, edges, wts, _bg);

    return true;
}

void computeWeights(const vector<point>& _nodes, const vector<ivec2>& _edges,
                    const vector<vector<float>>& _msures, Weight wt_type,
                    vector<float>& _weights)
{
    switch (wt_type)
    {
        case graphapp::BT3:
        {
            // weight is computed by inversion & normalization
            _weights = _msures[BT3];
            auto min_max =
                std::minmax_element(_weights.begin(), _weights.end());
            auto min = *min_max.first;
            auto max = *min_max.second;
            for (auto& v : _weights)
            {
                v = util::rescale(-v, -max, -min, 0.0f, 1.0f);
            }
            break;
        }
        default:
            _weights = _msures[BT3];
            break;
    }
}

bool makeBoostGraph(const vector<point>& _nodes, const vector<ivec2>& _edges,
                    const vector<float>& _weights, WeightedGraph& _g)
{
    for (auto i = 0; i < _nodes.size(); ++i)
        _g.add_vertex({i, _nodes[i]});

    for (auto i = 0; i < _edges.size(); ++i)
    {
        auto e = _edges[i];
        NodeHandle u = boost::vertex(e[0], _g);
        NodeHandle v = boost::vertex(e[1], _g);
        _g.add_edge(u, v, {i, _weights[i]});
    }

    // sanity check: is there index property by default?
    auto idx_map = boost::get(boost::vertex_index, _g);
    NodeIter v_it, v_end;
    boost::tie(v_it, v_end) = boost::vertices(_g);
    bool matched = true, indexable_by_outside_idx = true;
    for (auto i = 0; i < _nodes.size(); ++i)
    {
        auto v = boost::vertex(i, _g);
        auto idx = idx_map[v];
        if (i != idx)
            indexable_by_outside_idx = false;
        if (idx != _g[v].idx)
            matched = false;
    }
    if (!indexable_by_outside_idx || !matched)
    {
        if (!indexable_by_outside_idx)
        {
            cout << "WARNING: sanity checking failed. "
                 << "outside index doesn't correspond to 1st arg in "
                    "boost::vertex()!"
                 << endl;
        }
        if (!matched)
        {
            cout << "WARNING: sanity checking failed. "
                 << "index from index-map and that stored within vertex "
                    "doesn't match!"
                 << endl;
        }
        return false;
    }
    return true;
}

bool findRoot(const WeightedGraph& _g,
              const vector<vector<float>>& _edge_msures, NodeHandle& _root)
{
    // determine root by the vertex with largest radius
    //
    vector<float> radii_v(boost::num_vertices(_g), 0.0f);
    const auto& radii_e = _edge_msures[BT3];

    // compute radius for each vertex
    // boost::graph_traits<WeightedGraph>::edge_iterator e_it, e_end_it;
    EdgeIter e_it, e_end_it;
    for (boost::tie(e_it, e_end_it) = boost::edges(_g); e_it != e_end_it;
         ++e_it)
    {
        auto e_r = radii_e[_g[*e_it].idx];
        auto u = boost::source(*e_it, _g);
        auto v = boost::target(*e_it, _g);
        radii_v[_g[u].idx] += e_r;
        radii_v[_g[v].idx] += e_r;
    }

    // find the one with largest r
    NodeIter vit, vit_end;
    float max_r = 0.0f;
    for (boost::tie(vit, vit_end) = boost::vertices(_g); vit != vit_end; ++vit)
    {
        auto& r = radii_v[_g[*vit].idx];
        r /= (float)boost::degree(*vit, _g);
        if (r > max_r)
        {
            max_r = r;
            _root = *vit;
        }
    }

    return true;
}
void make_shortest_path_tree(const WeightedGraph& _g, const NodeHandle& _s,
                             vector<NodeHandle>& _parent);
void makeTreeFromGraph(const WeightedGraph& _g,
                       const vector<vector<float>>& _edge_msures,
                       TreeMethod _method, vector<NodeHandle>& _t)
{
    cout << "finding root ... " << endl;
    NodeHandle _root;
    if (!findRoot(_g, _edge_msures, _root))
    {
        cout << "Using default root (0th vertex)" << endl;
        _root = boost::vertex(0, _g);
    }
    else
    {
        cout << "Done: vertex " << _g[_root].idx << ", " << _g[_root].p
             << " is root." << endl;
    }

    switch (_method)
    {
        case graphapp::ShortestPath:
            make_shortest_path_tree(_g, _root, _t);
            break;
        default:
            make_shortest_path_tree(_g, _root, _t);
            break;
    }
}

void convertToEdges(const WeightedGraph& _g, const vector<NodeHandle>& _t,
                    vector<ivec2>& _edge_list)
{
    std::unordered_set<ivec2, ivec2Hash> edge_set;
    for (auto i = 0; i < boost::num_vertices(_g); ++i)
    {
        auto v_hdl = boost::vertex(i, _g);
        auto cur_idx = _g[v_hdl].idx;
        auto par_idx = _g[_t[i]].idx;
        if (par_idx != cur_idx)
        {
            edge_set.insert(util::makeEdge(cur_idx, par_idx));
        }
    }
    _edge_list.clear();
    _edge_list.insert(_edge_list.end(), edge_set.begin(), edge_set.end());
}

void make_shortest_path_tree(const WeightedGraph& _g, const NodeHandle& _s,
                             vector<NodeHandle>& _parent)
{
    try
    {
        _parent.resize(boost::num_vertices(_g));
        vector<float> dist(boost::num_vertices(_g));

        cout << "running boost dijkstra ... " << endl;
        // call boost shortest path tree.
        auto v_id_map = boost::get(boost::vertex_index, _g);
        auto e_wt_map = boost::get(&GraphEdge::w, _g);
        boost::dijkstra_shortest_paths(
            _g, _s,
            boost::weight_map(e_wt_map)
                .predecessor_map(boost::make_iterator_property_map(
                    _parent.begin(), v_id_map))
                .distance_map(
                    boost::make_iterator_property_map(dist.begin(), v_id_map)));
        cout << "Done: dijkstra completed." << endl;
    }
    catch (const boost::negative_edge& _ex)
    {
        cout << _ex.what() << endl;
    }
    /*catch ( const boost )*/
}

void exportTree(const WeightedGraph& _g, vector<NodeHandle>& _t,
                const std::string& _filename)
{
    vector<ivec2> edges;
    convertToEdges(_g, _t, edges);
    vector<trimesh::point> vts(boost::num_vertices(_g));
    for (auto i = 0; i < vts.size(); ++i)
    {
        auto p = _g[boost::vertex(i, _g)].p;
        vts[i].x = p[0];
        vts[i].y = p[1];
        vts[i].z = p[2];
    }
    auto err =
        voxelvoro::writeToPLY(_filename.c_str(), vts, edges, {}, {}, {}, {});
    if (err == voxelvoro::ExportErrCode::FAILURE)
    {
        throw std::runtime_error(std::string("Cannot write graph to file -> ") +
                                 _filename);
    }
    // ply::PLYWriter writer;
    // writer.write(_filename.c_str(), true, true, false, vts, edges, {});
}
} // namespace graphapp
