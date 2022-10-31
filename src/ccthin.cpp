#include "ccthin.h"

#include <cassert>
#include <iostream>
#include <list>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <vector>

using std::cout;
using std::endl;
using std::unordered_map;

namespace
{
// <out vert idx, rel order for sorting, poly dual or not>
typedef std::tuple<unsigned, float, bool> out_item;
void print_direct_g(const map<unsigned, list<out_item>>& _directed_graph);
void print_int_list(const vector<int>& _l);
} // namespace

CellComplexThinning::CellComplexThinning() { reset(); }

CellComplexThinning::~CellComplexThinning() {}

void CellComplexThinning::setup(cellcomplex* _cc) { m_cc = _cc; }

void CellComplexThinning::reset()
{
    m_cc = nullptr;

    /* information for pruning */
    // values used in pruning for HS elements
    for (auto& msre : m_measure)
    {
        msre.clear();
        msre.shrink_to_fit();
    }

    // prune specific ref count
    // valid after each prune
    for (auto rmlist : m_removed)
    {
        rmlist.clear();
        rmlist.shrink_to_fit();
    }
    m_ref_vert_per_prune.clear();
    m_ref_vert_per_prune.shrink_to_fit();
    m_ref_edge_per_prune.clear();
    m_ref_edge_per_prune.shrink_to_fit();
    // faces to remove flag
    m_to_remove_face.clear();
    m_to_remove_face.shrink_to_fit();
    m_face_area.clear();
    m_face_area.shrink_to_fit();
}

void CellComplexThinning::preprocess()
{
    ///*compute area of each face*/
    // float a_min = std::numeric_limits<float>::max();
    // float a_max = -1.0f;
    // m_face_area.resize( m_tri_faces.size() );
    // for ( size_t fi = 0; fi < m_tri_faces.size(); ++fi )
    //{
    //	const auto& f = m_tri_faces[ fi ];
    //	float area = trimesh::len(
    //		trimesh::trinorm( m_vts[ f[ 0 ] ], m_vts[ f[ 1 ] ], m_vts[ f[ 2
    //]
    //]
    //)
    //		);
    //	m_face_area[ fi ] = area;
    //	a_min = std::min( area, a_min );
    //	a_max = std::max( area, a_max );
    // }
    // cout << "skeleton: each face's area computed. [" << a_min << "," << a_max
    // << "]." << endl;

    /* prepare for pruning */

    // init remove tags
    m_removed[FACE].assign(m_cc->numFaces(), false);
    m_removed[EDGE].assign(m_cc->numEdges(), false);
    m_removed[VERTEX].assign(m_cc->numVts(), false);
    cout << "skeleton: removed-tags init.ed." << endl;

    // compute face area
    // m_face_area.reserve( m_cc->numFaces() );
    // vector<int> f_v_rep;
    // vector<uTriFace> tris;
    // for ( auto fi = 0; fi < m_cc->numFaces(); ++fi )
    //{
    //	m_cc->getFaceVRep( fi, f_v_rep );
    //	if ( f_v_rep.size() > 3 )
    //	{
    //		tris.clear();
    //		util::simpleTriangulate( f_v_rep, tris );
    //	}
    //	else
    //	{
    //		tris.push_back( uTriFace( f_v_rep.data() ) );
    //	}
    //	// sum up all triangles
    //	float area = 0.0f;
    //	for ( const auto& t : tris )
    //	{
    //		auto v1 = m_cc->getVert( t[ 1 ] ) - m_cc->getVert( t[ 0 ] );
    //		auto v2 = m_cc->getVert( t[ 2 ] ) - m_cc->getVert( t[ 0 ] );
    //		area += trimesh::len( v1.cross( v2 ) );
    //	}
    //	m_face_area.push_back( area );
    //}
}

const vector<bool>& CellComplexThinning::getEdgeRemovedTags() const
{
    return m_removed[EDGE];
}

const vector<bool>& CellComplexThinning::getFaceRemovedTags() const
{
    return m_removed[FACE];
}

void CellComplexThinning::getRemainedVts(vector<int>& _remain_vts_ids) const
{
    // the set of vertices of interest here are those vertices
    // referred to by at least one edge originally
    _remain_vts_ids.clear();
    for (size_t i = 0; i < m_cc->numVts(); ++i)
    {
        if (!m_removed[VERTEX][i])
        {
            _remain_vts_ids.push_back(i);
            /*if (m_ref_vert_per_prune[i] <= 0)
            {
                    cout << "Attention: vert "<<i<<" is not removed but has
            ref-cnt
            "<<m_ref_vert_per_prune[i]<<endl;
            }*/
        }
        else
        {
            if (m_ref_vert_per_prune[i] > 0)
            {
                cout << "Logic Error: vert " << i
                     << " is removed but still has ref-cnt "
                     << m_ref_vert_per_prune[i] << endl;
            }
        }
    }
}

void CellComplexThinning::getRemainedEdges(vector<int>& _edges_ids) const
{
    _edges_ids.clear();
    for (unsigned i = 0; i < m_cc->numEdges(); ++i)
    {
        // how about just returning remaining isolated edges,
        // i.e. those w/o being used by any face
        if (!m_removed[EDGE][i] && this->m_ref_edge_per_prune[i] == 0)
            _edges_ids.push_back(i);
    }
}

void CellComplexThinning::getRemainedFaces(vector<int>& _face_ids) const
{
    _face_ids.clear();
    for (unsigned fi = 0; fi < m_cc->numFaces(); ++fi)
    {
        if (!m_removed[FACE][fi])
            _face_ids.push_back(fi);
    }
}

void CellComplexThinning::assignElementValues(
    const vector<float>& _vert_measures, const vector<float>& _edge_measures,
    const vector<float>& _face_measures)
{
    m_measure[VERTEX] = _vert_measures;
    m_measure[EDGE] = _edge_measures;
    m_measure[FACE] = _face_measures;
    min_f_measure = numeric_limits<float>::max();
    max_f_measure = numeric_limits<float>::min();
    min_e_measure = min_f_measure;
    max_e_measure = max_f_measure;
}

void CellComplexThinning::setPolyFaceDegenerateThresh(float _t)
{
    m_face_degen_thresh = _t;
}

void CellComplexThinning::setComponentFaceNumberThresh(int _t)
{
    m_cmpnt_num_faces_thresh = _t;
}

void CellComplexThinning::prune(float _f_t, float _l_t,
                                bool _remove_small_components)
{
    auto sanity_check = [&]() {
        for (unsigned ei = 0; ei < m_cc->numEdges(); ++ei)
        {
            if (this->m_ref_edge_per_prune[ei] == 0 && !m_removed[0][ei])
            {
                const auto& e = m_cc->getEdge(ei);
                /*cout << "0-ref edge " << ei << e << " not removed! "
                        << ( part_of_orig_edge( ei ) ? "st edge." : "dual edge."
                   ) <<" "
                        << edges_measure[ei] <<". "
                        << "ref of vert at end: "
                        << m_ref_vert_per_prune[ e[ 0 ] ] << "," <<
                   m_ref_vert_per_prune[ e[ 1 ] ]
                        << endl;*/
                // break;
            }
        }
    };

    std::cout << "f_t, e_t: " << _f_t << ", " << _l_t << std::endl;

    set<unsigned> vts_to_debug;
    // vts_to_debug.insert(3807);
    queue<simple_pair> q;

    /*
    01. build reference count map
    02. push all simple pairs in q that are blow the thresholds
    03. iteratively retracting
    */
    /*01. reset remove tags & ref count for verts/edges */
    m_removed[EDGE].assign(m_cc->numEdges(), false);
    m_removed[FACE].assign(m_cc->numFaces(), false);
    m_removed[VERTEX].assign(m_cc->numVts(), false);
    m_ref_vert_per_prune = m_cc->refCntPerVert();
    m_ref_edge_per_prune = m_cc->refCntPerEdge();
    // for now no faces will be removed
    m_to_remove_face.assign(m_cc->numFaces(), false);

    /*02. init q: find all simple pairs whose simple cell has a value below
     * threshold */
    cout << "init.ing q ..." << endl;
    // push simple face-edge pairs to q
    for (unsigned ei = 0; ei < m_cc->numEdges(); ++ei)
    {
        if (m_ref_edge_per_prune[ei] == 1)
        {
            unsigned nb_f = m_cc->nbFaceofEdge(ei, 0);
            /*float val_bt2bt3 = edges_diff[ei];
            float val_bt2bt3rel = edges_reldiff[ei];*/
            if (face_edge_pair_below_threshold(nb_f, ei, _f_t, _l_t))
            {
                q.push(simple_pair(simple_pair::FE_PAIR, nb_f, ei));
            }
        }
    }
    // push simple edge-vertex pairs to q
    for (unsigned vi = 0; vi < m_cc->numVts(); ++vi)
    {
        if (m_ref_vert_per_prune[vi] == 1)
        {
            unsigned nb_e = m_cc->nbEdgeofVert(vi, 0);
            if (edge_vert_pair_below_threshold(nb_e, _l_t))
                q.push(simple_pair(simple_pair::EV_PAIR, nb_e, vi));
        }
    }
    cout << "after init, q size: " << q.size() << "" << endl;
    // cout << "prune: preparation done. " << endl;

    // vts_to_debug.insert(424);
    for (unsigned vi = 0; vi < m_cc->numVts(); ++vi) // debug
    {
        if (vts_to_debug.count(vi))
        {
            cout << "nbs of " << vi << ": ";
            print_int_list(m_cc->nbEdgesofVert(vi));
            cout << endl;
        }
    }

    prune_while_iteration(vts_to_debug, _f_t, _l_t, q);

    // removed and ref count should be consistent
    /*cout << "sanity check - after pruning, before non-dual edge removal ..."
    << endl; sanity_check();*/

    // removed and ref count should be consistent
    /*cout << "sanity check - after non-dual edge removal, before small
    components removal..." << endl; sanity_check();*/

    /* we may want to remove those small isolated components for better look */
    // if (_remove_small_components)
    //{
    //	cout << "removing *small* components..." << endl;

    //	// identify *small components* and remove all poly faces
    //	// and non-dual edges, and all relevant vertices
    //	mark_components( m_face_degen_thresh, m_cmpnt_num_faces_thresh );
    //	//cout << "small components marked!" << endl;

    //	// need to re-initialize q with new simple pairs
    //	for (unsigned ei = 0; ei < m_cc->numFaces(); ++ei)
    //	{
    //		if ( m_ref_edge_per_prune[ ei ] == 1 && !m_removed[ 0 ][ ei ] )
    //		{
    //			unsigned nb_f;
    //			for ( auto nb_f_i = 0; nb_f_i < m_cc->cntNbFacesofEdge(
    // ei
    //);
    //++nb_f_i )
    //			{
    //				nb_f = m_cc->nbFaceofEdge( ei, nb_f_i );
    //				if ( !m_removed[ FACE ][ nb_f ] )
    //				{
    //					break;
    //				}
    //			}
    //			if (face_edge_pair_below_threshold(nb_f, ei, _f_t,
    //_l_t))
    //			{
    //				q.push(simple_pair(simple_pair::FE_PAIR, nb_f,
    // ei));
    //			}
    //		}
    //	}
    //}
    // cout << "after marking small components, q size: " << q.size() << endl;

    // we may want to perform another pruning
    /*if ( !q.empty() )
    {
            prune_while_iteration(vts_to_debug, _f_t, _l_t, q);
    }*/

    // removed and ref count should be consistent
    /*cout << "sanity check - after small components removal..." << endl;
    sanity_check();*/

    /* post pruning sanity check */
    // removed and ref count should be consistent
    /*cout << "sanity check - after (clean-up) non-dual edge removal..." <<
    endl; sanity_check();*/

    // report what's left: any simple pair that should be removed however left?
    // first, edge-vert pair
    for (unsigned ei = 0; ei < m_cc->numEdges(); ++ei)
    {
        if (m_removed[0][ei])
            continue;

        const auto& e = m_cc->getEdge(ei);
        if (m_ref_edge_per_prune[ei] == 1 &&
            edge_vert_pair_below_threshold(ei, _l_t))
        {
            if (m_ref_vert_per_prune[e[0]] == 1)
            {
                cout << "edge-vert pair " << ei << "-" << e[0]
                     << "should be removed!" << endl;
            }
            if (m_ref_vert_per_prune[e[1]] == 1)
            {
                cout << "edge-vert pair " << ei << "-" << e[1]
                     << "should be removed!" << endl;
            }
        }
    }
    // second, examine face-edge pair
    vector<int> f;
    for (unsigned fi = 0; fi < m_cc->numFaces(); ++fi)
    {
        if (m_removed[1][fi])
            continue;

        m_cc->getFaceERep(fi, f);
        for (auto ei : f)
        {
            if (m_ref_edge_per_prune[ei] == 1 && !m_removed[EDGE][ei] &&
                !m_removed[FACE][fi] &&
                face_edge_pair_below_threshold(fi, ei, _f_t, _l_t))
            {
                cout
                    << "face-edge pair " << fi << "-" << ei
                    << "should have been removed!" << endl
                    << "(check simple pair removal logic, or the logic used to "
                       "perform this test.)"
                    << endl;
                // goto END_OF_POST_TEST;
            }
        }
    }
    // third, those faces marked as TO-REMOVE must have been removed
    for (unsigned fi = 0; fi < m_cc->numFaces(); ++fi)
    {
        if (m_to_remove_face[fi] && !m_removed[FACE][fi])
        {
            cout << "face " << fi << " (marked as to-remove) but not removed!"
                 << "value: " << m_measure[FACE][fi] << endl;
            m_cc->getFaceERep(fi, f);
            for (auto ei : f)
            {
                if (!m_removed[EDGE][ei])
                {
                    cout << "edge " << ei << " removed? " << m_removed[0][ei]
                         << ", "
                         << " ref count " << m_ref_edge_per_prune[ei] << ", "
                         << " remaining nb faces: ";
                    for (auto nb_f_i = 0; nb_f_i < m_cc->cntNbFacesofEdge(ei);
                         ++nb_f_i)
                    {
                        auto nb_f = m_cc->nbFaceofEdge(ei, nb_f_i);
                        if (!m_removed[FACE][nb_f])
                            cout << nb_f
                                 << (m_to_remove_face[nb_f] ? "(m)" : "")
                                 << ", ";
                    }
                    cout << endl;
                }
            }
            // goto END_OF_POST_TEST;
            // break;
        }
    }
END_OF_POST_TEST:;
}

cellcomplex CellComplexThinning::remainingCC()
{
    // remaining edges
    vector<int> edge_indices;
    vector<ivec2> edges_cc;
    getRemainedEdges(edge_indices);
    for (auto i : edge_indices)
    {
        edges_cc.push_back(m_cc->getEdge(i));
    }
    // remaining faces. will be broken into triangles.
    vector<int> face_indices;
    vector<uTriFace> tri_faces_cc;
    getRemainedFaces(face_indices);
    vector<int> f;
    vector<uTriFace> f_tris;
    for (auto i : face_indices)
    {
        m_cc->getFaceVRep(i, f);
        util::simpleTriangulate(f, f_tris);
        tri_faces_cc.insert(tri_faces_cc.end(), f_tris.begin(), f_tris.end());
    }
    // remaining vertices.
    vector<int> vts_remained_indices;
    vector<point> vts_cc;
    getRemainedVts(vts_remained_indices);
    for (auto i : vts_remained_indices)
    {
        vts_cc.push_back(m_cc->getVert(i));
    }

    cout << "# remaining vts ( >= # connected components): "
         << vts_remained_indices.size() << endl;
    cout << "# remaining edges: " << edge_indices.size() << endl;
    cout << "# remaining faces: " << face_indices.size() << endl;

    // compact remaining cc
    util::compactify(vts_remained_indices, edges_cc, tri_faces_cc);
    return cellcomplex(vts_cc, edges_cc, tri_faces_cc);
}

void CellComplexThinning::remainingCC(vector<point>& vts_cc,
                                      vector<ivec2>& edges_cc,
                                      vector<uTriFace>& _tri_faces,
                                      vector<int>* _from_fi_ptr) const
{
    // remaining edges
    vector<int> edge_indices;
    edges_cc.clear();
    getRemainedEdges(edge_indices);
    for (auto ei : edge_indices)
    {
        edges_cc.push_back(m_cc->getEdge(ei));
    }
    // remaining faces. will be broken into triangles.
    vector<int> face_indices;
    _tri_faces.clear();
    if (_from_fi_ptr)
        _from_fi_ptr->clear();
    getRemainedFaces(face_indices);
    vector<int> f;
    vector<uTriFace> f_tris;
    for (auto fi : face_indices)
    {
        m_cc->getFaceVRep(fi, f);
        util::simpleTriangulate(f, f_tris);
        _tri_faces.insert(_tri_faces.end(), f_tris.begin(), f_tris.end());
        if (_from_fi_ptr)
            _from_fi_ptr->insert(_from_fi_ptr->end(), f_tris.size(), fi);
    }
    // remaining vertices.
    vector<int> vts_remained_indices;
    vts_cc.clear();
    getRemainedVts(vts_remained_indices);
    for (auto i : vts_remained_indices)
    {
        vts_cc.push_back(m_cc->getVert(i));
    }

    cout << "# remaining vts after thinning ( >= # connected components): "
         << vts_remained_indices.size() << endl;
    cout << "# remaining edges: " << edge_indices.size() << endl;
    cout << "# remaining faces: " << face_indices.size() << endl;

    // lastly compact remaining cc
    util::compactify(vts_remained_indices, edges_cc, _tri_faces);
}

bool CellComplexThinning::edge_vert_pair_below_threshold(unsigned _ei,
                                                         float _t) const
{
    return m_measure[EDGE][_ei] < _t;
}

bool CellComplexThinning::face_edge_pair_below_threshold(unsigned _fi,
                                                         unsigned _ei,
                                                         float _f_t,
                                                         float _e_t) const
{
    return m_to_remove_face[_fi] ||
           m_measure[FACE][_fi] < _f_t /*&& m_measure[ EDGE ][ _ei ] < _e_t*/;
}

void CellComplexThinning::prune_while_iteration(
    const set<unsigned>& vts_to_debug, float _f_t, float _l_t,
    std::queue<simple_pair>& q)
{
    int ei_debug = -1; // debug
    while (!q.empty())
    {
        const auto spr = q.front();
        q.pop();

        // debug
        if (spr.type == simple_pair::EV_PAIR && vts_to_debug.count(spr.idx1))
        {
            cout << "simple pair popped: (" << spr.type << ", " << spr.idx0
                 << ", " << spr.idx1 << ")"
                 << "ref_vert: " << m_ref_vert_per_prune[spr.idx1] << endl;
        }
        if (spr.type == simple_pair::FE_PAIR && spr.idx1 == ei_debug)
        {
            cout << "got face-edge spr (" << spr.idx0 << "," << spr.idx1 << ")"
                 << " face removed?"
                 << (m_removed[FACE][spr.idx0] ? "true" : "false") << endl;
            print_edge(ei_debug);
        }

        // skip if this simple pair is removed already (the simple cell removed)
        if (m_removed[spr.type][spr.idx0])
        {
            assert(spr.type == 1 ? m_ref_edge_per_prune[spr.idx1] == 0
                                 : m_ref_vert_per_prune[spr.idx1] == 0);

            // debug
            if (spr.type == 0 && vts_to_debug.count(spr.idx1))
            {
                cout << "skip: simple cell removed." << endl;
            }

            continue;
        }

        // logic checkpoint
        if (!(spr.type == simple_pair::FE_PAIR &&
                  m_ref_edge_per_prune[spr.idx1] == 1 ||
              spr.type == simple_pair::EV_PAIR &&
                  m_ref_vert_per_prune[spr.idx1] == 1))
        {
            cout << "logic err: simple pair type " << spr.type;
            if (spr.type == simple_pair::FE_PAIR)
                cout << " ref_edge[" << spr.idx1
                     << "] = " << m_ref_edge_per_prune[spr.idx1];
            else
                cout << " ref_vert[" << spr.idx1
                     << "] = " << m_ref_vert_per_prune[spr.idx1];
            cout << endl;
            exit(1);
        }

        // still valid simple pair. remove it,
        // and modify ref count for the elements touched by the removal
        if (spr.type == simple_pair::EV_PAIR) // edge-vert pair
        {
            // set the simple pair: edge and vert removed
            m_removed[EDGE][spr.idx0] = true;
            m_removed[VERTEX][spr.idx1] = true;

            // modify ref count for vts of the removed edge
            const auto& e = m_cc->getEdge(spr.idx0);

            m_ref_vert_per_prune[e[0]]--;
            m_ref_vert_per_prune[e[1]]--;

            // debug
            if (vts_to_debug.count(spr.idx1))
            {
                cout << "after dec (due to edge-vert removal " << spr.idx0
                     << "-" << spr.idx1 << "), ref_vert[" << spr.idx1
                     << "]=" << m_ref_vert_per_prune[spr.idx1] << endl;
            }

            assert(m_ref_vert_per_prune[spr.idx1] == 0);

            // add the only new simple pair edge-vert to q
            // if a new simple pair arises
            unsigned other_end = e[0] != spr.idx1 ? e[0] : e[1];
            if (m_ref_vert_per_prune[other_end] == 1)
            {
                for (auto nb_e_i = 0;
                     nb_e_i < m_cc->cntNbEdgesofVert(other_end); ++nb_e_i)
                {
                    auto nb_e = m_cc->nbEdgeofVert(other_end, nb_e_i);
                    if (!m_removed[EDGE][nb_e])
                    {
                        if (edge_vert_pair_below_threshold(nb_e, _l_t))
                        {
                            q.push(simple_pair(simple_pair::EV_PAIR, nb_e,
                                               other_end));

                            if (vts_to_debug.count(other_end)) // debug
                            {
                                cout << "simple pair added: (" << 0 << ", "
                                     << nb_e << ", " << other_end << ")"
                                     << endl;
                            }
                        }
                        break;
                    }
                }
            }
        }
        else // face-edge pair
        {
            m_removed[FACE][spr.idx0] = true;
            m_removed[EDGE][spr.idx1] = true;

            /*if (spr.idx0 == 176 && spr.idx1 == 1388)
            {
            cout << spr.idx0<<"removed? "<<removed[1][spr.idx0]<<" "
            <<spr.idx1<<"removed? "<<removed[0][spr.idx1]<<endl;
            }*/

            m_ref_edge_per_prune[spr.idx1]--;
            assert(m_ref_edge_per_prune[spr.idx1] == 0);

            if (spr.idx1 == ei_debug)
            {
                cout << "after face-edge spr (" << spr.idx0 << "," << spr.idx1
                     << ") removed" << endl;
                print_edge(ei_debug);
            }

            // modify ref count for the un-removed edges of the removed face
            vector<int> f;
            m_cc->getFaceERep(spr.idx0, f);
            for (auto ei : f)
            {
                if (!m_removed[EDGE][ei])
                {
                    m_ref_edge_per_prune[ei]--;

                    if (spr.idx1 == ei_debug)
                    {
                        cout << "after modifying other edges of the face"
                             << endl;
                        print_edge(ei_debug);
                    }

                    // any new face-edge simple pair?
                    if (m_ref_edge_per_prune[ei] == 1)
                    {
                        // find the unique face
                        for (auto nb_f_i = 0;
                             nb_f_i < m_cc->cntNbFacesofEdge(ei); ++nb_f_i)
                        {
                            auto nb_f = m_cc->nbFaceofEdge(ei, nb_f_i);
                            if (!m_removed[FACE][nb_f])
                            {
                                // found the face. add to q the pair if criteria
                                // met
                                if (face_edge_pair_below_threshold(nb_f, ei,
                                                                   _f_t, _l_t))
                                {
                                    q.push(simple_pair(simple_pair::FE_PAIR,
                                                       nb_f, ei));
                                }
                                // if ( polyFace_bt2bt3[*nb_f] <= _bt2bt3_t ||
                                //	polyFace_bt2bt3rel[*nb_f] <=
                                //_bt2bt3rel_t )
                                //{
                                //	if (part_of_orig_edge(*eit))
                                //	{
                                //		if (edges_diff[*eit] <=
                                //_bt2bt3_t
                                //|| 			edges_diff[*eit] <=
                                //_bt2bt3rel_t)
                                //		{
                                //			q.push(simple_pair(1,
                                //*nb_f, *eit));
                                //		}
                                //	}
                                //	else
                                //	{
                                //		if (edges_diff[*eit] <=
                                //_bt1bt2_t
                                //|| 			edges_diff[*eit] <=
                                //_bt1bt2rel_t)
                                //		{
                                //			q.push(simple_pair(1,
                                //*nb_f, *eit));
                                //		}
                                //	}
                                // }
                                break;
                            }
                        }
                    }
                }
            }

            // modify ref count for vts of the edge of the removed edge
            auto e = m_cc->getEdge(spr.idx1);
            for (unsigned i = 0; i < 2; ++i)
            {
                unsigned vi = e[(int)i];
                m_ref_vert_per_prune[vi]--;
                if (vts_to_debug.count(vi)) // debug
                {
                    cout << "after dec (due to face-edge removal. type "
                         << spr.type << ", " << spr.idx0 << "-" << spr.idx1
                         << ", ref_vert[" << vi
                         << "]=" << m_ref_vert_per_prune[vi] << endl;
                }

                // any new edge-vert simple pair?
                if (m_ref_vert_per_prune[vi] == 1)
                {
                    // locate the unique nb edge
                    for (auto ith_nb_e = 0;
                         ith_nb_e < m_cc->cntNbEdgesofVert(vi); ++ith_nb_e)
                    {
                        auto nb_ei = m_cc->nbEdgeofVert(vi, ith_nb_e);
                        if (!m_removed[EDGE][nb_ei])
                        {
                            // edge located. add to q the pair if criteria met.
                            if (edge_vert_pair_below_threshold(nb_ei, _l_t))
                            {
                                q.push(simple_pair(simple_pair::EV_PAIR, nb_ei,
                                                   vi));

                                if (vts_to_debug.count(vi)) // debug
                                {
                                    cout << "simple pair added: (" << 0 << ", "
                                         << nb_ei << ", " << vi << ")" << endl;
                                }
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
    cout << "prune: iterative removal done. " << endl;
}

void CellComplexThinning::mark_components(float _cmpnt_geomSize_t,
                                          int _num_face_t)
{
    vector<bool> edge_visited(m_cc->numEdges(), false);

    /*
    while there is unvisited edge:
            1. find a cmpnt C; mark edges as visited in the process;
            2. if |C| has too many faces, mark all its faces as NOT-TO-REMOVE;
            3. else, mark all faces as REMOVE, if C is *small* enough;
    */
    queue<unsigned> q_faces;  // queue of faces' id.s
    set<unsigned> face_cmpnt; // store faces of the cur component uniquely
    vector<point> vts_pos_cmpnt;
    vector<int> f_e_rep;
    int cnt_total_cmpnts = 0;
    int cnt_small_cmpnts = 0;
    for (unsigned ei = 0; ei < m_cc->numEdges(); ++ei)
    {
        if (m_ref_edge_per_prune[ei] != 1)
            continue;

        if (edge_visited[ei])
            continue;

        edge_visited[ei] = true;

        // start growing from this edge:
        // initialize queue with this only face
        int seed_fi = -1;
        for (auto ith_nb_f = 0; ith_nb_f < m_cc->cntNbFacesofEdge(ei); ith_nb_f)
        {
            auto nb_fi = m_cc->nbFaceofEdge(ei, ith_nb_f);
            if (!m_removed[FACE][nb_fi])
            {
                seed_fi = nb_fi;
                break;
            }
        }
        assert(seed_fi >= 0);
        face_cmpnt.clear();
        face_cmpnt.insert(seed_fi);
        q_faces.push(seed_fi);

        // start growing from this face:
        f_e_rep.clear();
        while (!q_faces.empty())
        {
            auto fi = q_faces.front();
            q_faces.pop();

            m_cc->getFaceERep(fi, f_e_rep);

            // find each nb face of this face f (a list of edge indices)
            for (unsigned i = 0; i < 3 /*f.size()*/; ++i)
            {
                auto ei = f_e_rep[i];
                if (edge_visited[ei])
                    continue;

                edge_visited[ei] = true;
                // only inspect edges that are still used by other face(s)
                // if ( _ref_edge[f[i]] > 1 )
                {
                    for (auto ith_nb_f = 0;
                         ith_nb_f < m_cc->cntNbFacesofEdge(ei); ith_nb_f)
                    {
                        auto nb_fi = m_cc->nbFaceofEdge(ei, ith_nb_f);
                        // skip nb face that's already removed
                        if (m_removed[FACE][nb_fi] || nb_fi == fi)
                            continue;

                        // only push those nb faces that are not yet in the
                        // current component to the queue
                        auto pre_size = face_cmpnt.size();
                        face_cmpnt.insert(nb_fi);
                        q_faces.push(nb_fi);
                    }
                }
            } // for (): done grabbing nb faces of f to component and queue
        }     // while(): done for cur component

        cnt_total_cmpnts++;

        // now component is ready. let's find out if it can be pruned.
        // is it too big?
        if (face_cmpnt.size() > _num_face_t)
        {
            // component is too big (too many faces, can't be an isolated small
            // component, right?) don't remove it.
            continue;
        }
        else
        {
            // test if its span in the space is small enough (area)
            float sum_area = .0f;
            for (auto fi_iter = face_cmpnt.begin(); fi_iter != face_cmpnt.end();
                 ++fi_iter)
            {
                sum_area += m_face_area[*fi_iter];
            }
            bool is_geom_small = sum_area <= _cmpnt_geomSize_t;
            if (is_geom_small)
            {
                cnt_small_cmpnts++;
                // this component is geometrically small.
                // mark all its faces as TO-REMOVE!
                for (auto fi_iter = face_cmpnt.begin();
                     fi_iter != face_cmpnt.end(); ++fi_iter)
                {
                    m_to_remove_face[*fi_iter] = true;
                }
            }
        } // if... else... : done testing whether to keep this component or not.
    }     // for each edge: done looping thru each edge to locate disconnected
          // components

    cout << cnt_total_cmpnts << " total components found. " << endl;
    cout << cnt_small_cmpnts
         << " small components found (marked as to-remove). " << endl;

} // mark_components()

void CellComplexThinning::print_edge(int _ei) const
{
    const auto& e = m_cc->getEdge(_ei);
    cout << "debug ---" << endl;
    cout << "edge: " << _ei << ", " << e << ". "
         << "removed? " << (m_removed[VERTEX][_ei] ? "true" : "false") << ". "
         << "#ref-by-face " << m_ref_edge_per_prune[_ei] << ". " << endl;
    cout << "--- debug " << endl;
}

// functions local to this translation unit
namespace
{
//
// debug & test routines
//
void print_direct_g(const map<unsigned, list<out_item>>& _directed_graph)
{
    for (auto it = _directed_graph.begin(); it != _directed_graph.end(); ++it)
    {
        cout << it->first << ": ";
        const auto& out_nodes = it->second;
        for (auto nit = out_nodes.begin(); nit != out_nodes.end(); ++nit)
        {
            cout << "(" << std::get<0>(*nit) << "," << std::get<1>(*nit) << ","
                 << std::get<2>(*nit) << ") -> ";
        }
        cout << '\b' << '\b' << '\b';
        cout << endl;
    }
}

void print_int_list(const vector<int>& _l)
{
    for (auto it = _l.begin(); it != _l.end(); ++it)
        cout << *it << " ";
}
} // namespace