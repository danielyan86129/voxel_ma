#ifndef EDGE_COLLAPSE_H
#define EDGE_COLLAPSE_H

//////////////////////////////////////////////////////////////////////////
// This file contains definition for edge collapse operator(s)
//////////////////////////////////////////////////////////////////////////

#include "commondefs.h"
#include <memory>
#include <vector>

namespace voxelvoro
{
using std::shared_ptr;
using std::vector;

class EdgeCollapser
{
public:
    EdgeCollapser() {}
    ~EdgeCollapser() {}

    /// implement these in subclasses
    // do whatever preprocessing here
    virtual void preProcess() = 0;
    // collapse all edges that are collapsible
    virtual void collapse() = 0;

protected:
    // point to a vector of vts' coords
    shared_ptr<vector<point>> m_vts_ptr;
    // point to a vector of tri faces
    shared_ptr<vector<uTriFace>> m_tris_ptr;
    // point to a vector of tri edges
    shared_ptr<vector<ivec2>> m_edges;
};

class TopoPreservEdgeCollapser : public EdgeCollapser
{
public:
    //
    // @construct this topo-preserving edge collapser from vertices and faces
    // of a triangle mesh (could be non-manifold)
    //
    TopoPreservEdgeCollapser(const shared_ptr<vector<point>>& _vts_ptr,
                             const shared_ptr<vector<uTriFace>>& _tris_ptr);
    ~TopoPreservEdgeCollapser();

    void preProcess();
    void collapse();

protected:
    /*helpers*/
    void compute_order();

    // the order (the topo complexity) of each vertex
    vector<int> m_vert_ord;
    // the order (the topo complexity) of each edge
    vector<int> m_edge_ord;
};
} // namespace voxelvoro

#endif // EDGE_COLLAPSE_H