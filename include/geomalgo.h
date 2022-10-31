#ifndef GEOM_ALGO_H
#define GEOM_ALGO_H

#include "commondefs.h"
#include <limits>
#include <vector>

namespace util
{
using std::vector;

// scale _v \in[min, max] to [0, 1]
template <typename T>
T unitize(T _v, T _min, T _max)
{
    auto range = std::max(0.0000000000001, double(_max - _min));
    return (T)((double(_v) - _min) / range);
}
// scale v \in [0,1] into new range [newmin, newmax]
template <typename T>
T rescale(T v, T new_min, T new_max)
{
    return (T)(double(v) * (new_max - new_min) + new_min);
}
// scale v from [oldmin, oldmax] to [newmin, newmax]
template <typename T>
T rescale(T v, T old_min, T old_max, T new_min, T new_max)
{
    return (T)(double(unitize(v, old_min, old_max)) * (new_max - new_min) +
               new_min);
}

// return the position of the highest bit that's set
inline int highestBit(int _n)
{
    if (!_n)
        return 0;
    int ret = 1;
    while (_n >>= 1)
        ret++;
    return ret - 1;
}

// log2 of an integer
inline int log2Int(int _n)
{
    auto hib = highestBit(_n);
    return hib + ((1 << hib) < _n ? 1 : 0);
}

/*
** utility: for floats
*/

//
// are two numbers "equal" up to a small difference?
//
template <typename T>
bool is_equal(T _a, T _b, T _eps)
{
    return std::abs(_a - _b) <= _eps;
}

//
// triangulate a polygon (made up of int indices) in the most simple way
// Note:
// - only *correct* when it's planar-convex
// - _poly's int indices must be in some order (CW or CCW)
//
void simpleTriangulate(const vector<int>& _poly, vector<uTriFace>& _tris);
//
// triangulate a polygon (made up of vertex coordinates)
// Note:
// - only *correct* when it's planar-convex
// - _poly's vertices must be in some order (CW or CCW)
//
void nonDegenTriangulate(const vector<point>& _poly, vector<uTriFace>& _tris);
// make edge from given 2 index numbers
inline ivec2 makeEdge(int _i, int _j)
{
    return _i > _j ? ivec2(_j, _i) : ivec2(_i, _j);
}

// make edge from a given face
inline void makeEdgesFromTri(const uTriFace& _triFace, ivec2* _edges)
{
    for (int i = 0; i < 3; ++i)
    {
        _edges[i] = makeEdge(_triFace[i], _triFace[(i + 1) % 3]);
    }
}

// return the vertex id opposite to the given edge in the given face
inline int oppoVertex(const uTriFace& _triFace, const ivec2& _e)
{
    if (makeEdge(_triFace[0], _triFace[1]) == _e)
        return _triFace[2];
    if (makeEdge(_triFace[0], _triFace[2]) == _e)
        return _triFace[1];
    return _triFace[0];
}
// return the edge opposite to the given vertex id in the given face
inline ivec2 oppoEdge(const uTriFace& _triFace, int _vi)
{
    if (_triFace[0] == _vi)
        return makeEdge(_triFace[1], _triFace[2]);
    if (_triFace[1] == _vi)
        return makeEdge(_triFace[0], _triFace[2]);
    return makeEdge(_triFace[0], _triFace[1]);
}

// test whether two given vertices are close
inline bool closeVertices(const point& _v1, const point& _v2, float _eps)
{
    return trimesh::len2(_v1 - _v2) <= _eps * _eps;
}

// replace the vertex id i with the given id j in face f
// if i is not found in f, then nothing changes
inline void replaceVinF(uTriFace& _f, int _i, int _j)
{
    _f[0] == _i ? _f[0] = _j
                : _f[1] == _i ? _f[1] = _j : _f[2] == _i ? _f[2] = _j : 0;
}

// trace face boundary (a seq of vertices) into a seq of edges
inline void traceFace(const vector<int>& _f_of_vts, vector<ivec2>& _f_of_edges)
{
    _f_of_edges.resize(_f_of_vts.size());
    for (size_t i = 0; i < _f_of_vts.size(); ++i)
    {
        int u = _f_of_vts[i];
        int v = _f_of_vts[(i + 1) % _f_of_vts.size()];
        _f_of_edges[i] = util::makeEdge(u, v);
    }
}

//
// whether the input polygon is *almost* degenerate to a line
// @param the extent of *almost* is determined by eps
//
bool is_degenerate(const vector<point> _poly, float _eps);

// find the num of connected components
int findNumConnComponents(const vector<point>& _vts,
                          const vector<ivec2>& _edges);
// user provides vert-vert adj list, a pre-allocated list for visited flag
int findNumConnComponents(const vector<vector<int>>& _v_v_adj_tbl,
                          const vector<point>& _vts, vector<bool>& _visited);
int findNumConnComponents( // TODO: get rid of the parameter: _n_vts
    const vector<vector<int>>& _v_v_adj_tbl, int _n_vts,
    vector<bool>& _visited);

// compute the euler characteristic number for the given triangle faces
// representing a mesh
void computeEulerChar(const vector<uTriFace>& _tris, eulerchar& _euler_struct);
void computeEulerChar(const vector<point>& _vts, const vector<ivec2>& _edges,
                      const vector<uTriFace>& _tris, eulerchar& _euler_struct);

// given a vertex-id list (a sub-sequence of some full vertex list),
// compact-ify the list and rename the edge list and face list accordingly
// to refer to the compact v-list, essentially the following effect:
// e.g.
// v-list {4, 5, 3, 9} -> {0, 1, 2, 3} conceptually
// e-list {<3, 5>, <3, 9>, <5, 9>} -> {<2, 1>, <2, 3>, <1, 3>} actual change
// f-list {<3, 5, 9>} -> {<2, 1, 3>} actual change
void compactify(const vector<int>& _v_l, vector<ivec2>& _e_l,
                vector<uTriFace>& _f_l);

//
// Return a RGB colour value (hot-cold scheme) given a scalar v in the range
// [vmin,vmax] In this case each colour component ranges from 0 (no
// contribution) to 1 (fully saturated), modifications for other ranges is
// trivial. The colour is clipped at the end of the scales if v is outside the
// range [vmin,vmax]
TriColor GetColour(float v, float vmin, float vmax);

//
// return whether two axis-aligned boxes (origin + dimensions) intersect
bool intersect(const ivec3& _o1, const ivec3& _dim1, const ivec3& _o2,
               const ivec3& _dim2);
} // namespace util

#endif // GEOM_ALGO_H