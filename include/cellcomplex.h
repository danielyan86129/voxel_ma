#pragma once
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

using std::ifstream;
using std::iostream;
using std::map;
using std::shared_ptr;
using std::string;
using std::vector;

#include "commondefs.h"
#include "geomalgo.h"
#include "spaceinfo.h"

/************************************************************************
** a bare - bone struct containing points, edges, and polygons
** fill this struct up with info from your own algorithm
************************************************************************/
struct cellcomplex
{
    /*
    ** constructors & destructors
    */
    cellcomplex();

    // construct from given v, e, and triangles.
    // by default this cc will be finalized. set _perform_finalize = false to
    // disable this feature. This minimizes memory footprint.
    cellcomplex(const vector<point>& _vts, const vector<ivec2>& _edges,
                const vector<uTriFace>& _faces, bool _perform_finalize = true);

    ~cellcomplex();

    /*
    ** interfaces for making query
    */
    inline bool isFinalized() const;
    // TODO: maybe we need to be able to know manifold-ness
    // bool is2Manifold();

    //
    // info will be computed if not present upon request
    // return true if computation takes place, false if no computation needs to
    // be done
    bool needVVAdjacency();

    //
    // collect edges from faces internally
    // return true if all edges are just computed, false if no computation needs
    // to be done
    bool needEdges();

    //
    // reset states
    void invalidateStates();

    //
    // compute connected components based on face adjacency, i.e.
    // only consider two faces sharing an edge as connected
    // @param _subset_faces: set to nullptr if want to compute for all faces
    int compute_conn_cmpnts(const vector<int>* _subset_faces = nullptr) const;
    //
    // compute the num of connected components based on vertex adjacenty, i.e.
    // two vertices belong to same conn comp as long as connected by an edge
    // @param _subset_faces: set to nullptr if want to compute for all faces
    int compute_conn_cmpnts(const vector<int>* _subset_faces,
                            int overload_func_tag) const;
    //
    // pre-allocate space for vts, edges, faces
    //
    inline void allocVts(int _n)
    {
        m_vts.reserve(_n);
        m_nb_edges_for_v.reserve(_n);
    }
    inline void allocEdges(int _n)
    {
        m_edges.reserve(_n);
        m_nb_faces_for_e.reserve(_n);
    }
    inline void allocFaces(int _n)
    {
        m_faceVtsN.reserve(_n);
        m_faceIndex.reserve(_n);
    }
    //
    // collect unused space
    //
    inline void collectSpace()
    {
        m_vts.shrink_to_fit();
        m_nb_edges_for_v.shrink_to_fit();
        m_edges.shrink_to_fit();
        m_nb_faces_for_e.shrink_to_fit();
        m_faceVtsN.shrink_to_fit();
        m_faceIndex.shrink_to_fit();
    }
    //
    // finalize this cell complex: collecting all 1-cell (edges) from 2-cells
    // (faces), counting references, etc.
    //
    void finalize();
    //
    // compute euler characteristics for this complex
    void eulerChar(eulerchar& _ec) const;
    void eulerChar(eulerchar& _ec);
    //
    // returns primitive (vert, edge, or face) by its index
    //
    inline const point& getVert(size_t _vi) const { return m_vts[_vi]; }
    inline const ivec2& getEdge(size_t _ei) const { return m_edges[_ei]; }
    //
    // convert list-of-vertex-indices to a list-of-vertex-coords
    //
    inline void getVertCoords(const vector<int>& _vts_indices,
                              vector<point>& _vts_coords) const
    {
        _vts_coords.clear();
        _vts_coords.reserve(_vts_indices.size());
        for (auto vi : _vts_indices)
        {
            _vts_coords.push_back(m_vts[vi]);
        }
    }
    //
    // returns the list-of-vertex-indices of a face
    //
    inline void getFaceVRep(int _fi, vector<int>& _f_of_vts) const
    {
        const auto& n = m_faceVtsN[_fi];
        const auto& i = m_faceIndex[_fi];
        const auto& face_list = m_faces_v_rep.find(n)->second;
        _f_of_vts.assign(face_list.begin() + i, face_list.begin() + i + n);
    }
    //
    // returns the list-of-edge indices of a face
    //
    inline void getFaceERep(int _fi, vector<int>& _f_of_edges) const
    {
        const auto& n = m_faceVtsN[_fi];
        const auto& i = m_faceIndex[_fi];
        const auto& face_list = m_faces_e_rep.find(n)->second;
        _f_of_edges.assign(face_list.begin() + i, face_list.begin() + i + n);
    }
    //
    // returns the list-of-edges of a face
    //
    inline void getFaceERep(int _fi, vector<ivec2>& _f_of_edges) const
    {
        vector<int> f_of_vts;
        getFaceVRep(_fi, f_of_vts);
        util::traceFace(f_of_vts, _f_of_edges);
    }
    //
    // return the ref to the list of vertices
    //
    inline const vector<point>& getAllVts() const { return m_vts; }
    //
    // return ref count (by faces) for the given edge id
    //
    inline int refCount(size_t _ei) const { return m_edge_ref_cnt[_ei]; }
    //
    // num of faces
    //
    inline size_t numFaces() const { return m_faceIndex.size(); }
    //
    // num of vertices
    //
    inline size_t numVts() const { return m_vts.size(); }
    //
    // num of edges
    //
    inline size_t numEdges() const { return m_edges.size(); }
    //
    // add a face (a list of vertex indices. triangle is a special case.)
    //
    void appendFace(const vector<int>& _f_of_vts);
    void appendFace(const uTriFace& _f_tri);
    //
    // add a face, both its vertex rep. and edge rep.
    // Note: caller should make sure the vert-rep and edge-rep
    // must correspond to the same face
    //
    void appendFace(const vector<int>& _f_of_vts,
                    const vector<int>& _f_of_sides);
    //
    // add an edge (could be referenced by a face as part of its boundary, or an
    // isolated edge)
    //
    inline void appendEdge(const ivec2& _e) { m_edges.push_back(_e); }
    //
    // add a vertex
    //
    inline void appendVert(const point& _v) { m_vts.push_back(_v); }
    //
    // create the globally unique infinite element (vertex/edge)
    // and return the element index (within vertex/edge list)
    //
    size_t createInfiniteVertex();
    size_t createInfiniteEdge();
    size_t createInfiniteVertex(const point& _p);
    //
    // clears storage
    //
    void clear();
    void clearFaceList();
    void clearEdgeList();
    void clearVertList();
    /// various getters
    //
    // return nb vts of given vertex
    //
    inline void nbVts(unsigned _vi, vector<int>& _nbs) const;
    //
    // return the cnt of neighbor elements (nb edge or face)
    //
    inline int cntNbEdgesofVert(unsigned _vi) const;
    inline int cntNbFacesofEdge(unsigned _ei) const;
    //
    // return the i-th neighbor edge/face of the vertex/edge
    //
    inline int nbEdgeofVert(unsigned _vi, int _ith_nb_e) const;
    inline int nbFaceofEdge(unsigned _ei, int _ith_nb_f) const;
    //
    // return a list of nb edges/faces indices
    //
    inline vector<int> nbEdgesofVert(unsigned _vi) const;
    inline vector<int> nbFacesofEdge(unsigned _ei) const;
    //
    // return a list of ref-cnt for vert/edge by nb edges/faces
    //
    vector<int> refCntPerVert() const;
    vector<int> refCntPerEdge() const;

private:
    /*
    ** private helpers
    */
    /// adjacency management related helpers.
    /// Implementation details likely to change depending on the actual layout
    // initialize adjacency info and get ready for building the actual adjacency
    void init_V_E_adjacency();
    void init_E_F_adjacency();
    void init_adjacency();
    // add a neighbor edge to the vert's adj. info list
    void add_nb_edge(int _vi, int _ei);
    // add a neighbor face to the edge's adj. info list
    void add_nb_face(int _ei, int _nb_fi);
    // clear adj storage
    void clear_adjacency();

public:
    /*
    ** public data members
    */
    // 0-cells
    vector<point> m_vts;
    // 1-cells
    vector<ivec2> m_edges;
    // 1-cells ref-count (by 2-cells)
    vector<int> m_edge_ref_cnt;
    // polygon face info
    // How this works: all the faces with k vts (k = m_faceVtsN[i])
    // are held in m_faces_by_sides[k], i.e. a vector of vts indices
    // then m_faceIndex further determines the offset within m_faces_by_sides[k]
    map<int, vector<int>> m_faces_v_rep;
    map<int, vector<int>> m_faces_e_rep;
    vector<int> m_faceIndex;
    vector<int> m_faceVtsN;

private:
    /*
    ** private data members
    */
    /** states **/
    // is cc finalized
    bool m_is_finalized;
    // does this complex have infinite elements?
    bool m_has_infinite_v;
    bool m_has_infinite_e;
    // indices referring to the infinite elements
    size_t m_inf_v_id;
    size_t m_inf_e_id;
    // adjacency info: nb edges for vert
    vector<vector<int>> m_nb_edges_for_v;
    // adjacency info: nb faces for edge
    vector<vector<int>> m_nb_faces_for_e;
};

#include "cellcomplex_imp.h"