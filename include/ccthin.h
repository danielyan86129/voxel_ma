#ifndef HYBRID_SKELETON_H
#define HYBRID_SKELETON_H

#include "cellcomplex.h"
#include "commondefs.h"
#include <map>
#include <queue>
#include <set>
#include <trimesh/XForm.h> // trimesh's transformation matrix

using namespace std;

class CellComplexThinning
{
public:
    enum celltype
    {
        FACE = 1,
        EDGE = 0,
        VERTEX = 2
    };
    struct simple_pair
    {
        enum spairtype
        {
            FE_PAIR = 1,
            EV_PAIR = 0
        };
        simple_pair(spairtype _t, unsigned _idx0, unsigned _idx1)
        {
            type = _t;
            idx0 = _idx0;
            idx1 = _idx1;
        }
        simple_pair(const simple_pair& _spr)
        {
            type = _spr.type;
            idx0 = _spr.idx0;
            idx1 = _spr.idx1;
        }
        simple_pair& operator=(const simple_pair& _spr)
        {
            type = _spr.type;
            idx0 = _spr.idx0;
            idx1 = _spr.idx1;

            return *this;
        }
        // type = 1: face-edge pair, idx0 is face id, idx1 is edge id
        // type = 0: edge-vert pair, idx0 is edge id, idx1 is vert id
        spairtype type;
        unsigned idx0;
        unsigned idx1;
    };

public:
    CellComplexThinning();
    ~CellComplexThinning();

    //
    // initialize h.s. states
    //
    void setup(cellcomplex* _cc);
    //
    // reset this thinning operator
    //
    void reset();

    void preprocess();

    const vector<bool>& getEdgeRemovedTags() const;
    const vector<bool>& getFaceRemovedTags() const;

    //
    // those not removed vts (0-cell)
    //
    void getRemainedVts(vector<int>& _vts_ids) const;
    //
    // those isolated and not removed edges (1-cell)
    //
    void getRemainedEdges(vector<int>& _edges_ids) const;
    //
    // those not removed faces (2-cell)
    //
    void getRemainedFaces(vector<int>& _face_ids) const;

    //
    // assign values to elements: faces- edge- and vertex-cells
    // thinning will be performed based on these values
    //
    void assignElementValues(const vector<float>& _vert_measures,
                             const vector<float>& _edge_measures,
                             const vector<float>& _face_measures);

    //
    // set the threshold for detecting degenerate poly faces
    //
    void setPolyFaceDegenerateThresh(float _t);
    void setComponentFaceNumberThresh(int _t);

    //
    // perform thinning on c.c.
    //
    void prune(float _f_t, float _l_t, bool _remove_small_components);

    //
    // return the remaining part of cc after thinning as a cc.
    // faces will be triangulated, and the face index for each triangle is
    // recorded
    cellcomplex remainingCC();
    void remainingCC(vector<point>& vts_cc, vector<ivec2>& edges_cc,
                     vector<uTriFace>& tri_faces_cc,
                     vector<int>* _from_fi = nullptr) const;

private:
    //
    // if the edge-vert pair has a value below threshold
    //
    bool edge_vert_pair_below_threshold(unsigned _ei, float _t) const;
    //
    // if the face-edge pair has a value below threshold
    //
    bool face_edge_pair_below_threshold(unsigned _fi, unsigned _ei, float _f_t,
                                        float _e_t) const;
    //
    //
    //
    void prune_while_iteration(const set<unsigned>& vts_to_debug, float _f_t,
                               float _l_t, std::queue<simple_pair>& q);
    void mark_components(float _cmpnt_geomSize_t, int _num_face_t);

    //
    // debug routines
    //
    void print_edge(int _ei) const;

private:
    // basic information H.S. needs from the outside
    cellcomplex* m_cc;

    /* information for pruning */
    // threshold under which faces are considered as degenerate
    float m_face_degen_thresh;
    // threshold above which a component will be considered as *too big*
    int m_cmpnt_num_faces_thresh;
    // values used in pruning for HS elements
    vector<float> m_measure[3];
    float min_f_measure, max_f_measure;
    float min_e_measure, max_e_measure;
    // remove tag for face [1], edge [0], and vertex [2]
    vector<bool> m_removed[3];
    // prune-specific ref count
    // valid after each thinning
    vector<int> m_ref_vert_per_prune;
    vector<int> m_ref_edge_per_prune;
    // prune-specific parent element id
    // the id refers to the unique parent
    // TODO: determine whether this is needed
    vector<int> m_unique_par_for_v_per_prune;
    vector<int> m_unique_par_for_e_per_prune;
    // area for each face
    vector<float> m_face_area;

    // remove flag for faces.
    // could be modified by small-component-marking process
    vector<bool> m_to_remove_face;
};
#endif