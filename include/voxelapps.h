#pragma once
#include "commondefs.h"
#include "voroinfo.h"
#include <vector>

using std::vector;

namespace voxelvoro
{
namespace apps
{
/************************
** helpers
*************************/
//
//

//
// given vts, edges, and faces that define a 2-cell complex
// return a list of label (for vertex), each corresponding to a segment (broken
// at junction usually)
void segment(const vector<point>& _skel_vts, const vector<ivec2>& _skel_edges,
             const vector<uTriFace>& _skel_faces, vector<int>& _skel_vts_label);

//
// the measure type on MC:
enum MCMeasure
{
    BT3,
    BT2,
    RIDGE,
    LENGTH,
    SegmentLabel
};

//
// Given the medial-curve (MC) structure and the scalar field on it
// assign a proper value from the field to each site
void assignScalarToSites(const VoroInfo& _voro, const vector<point>& _mc_vts,
                         const vector<float>& _mc_msure,
                         const vector<int>& _mc_order, MCMeasure _msure_type,
                         float _smooth_ratio,
                         const CellComplexThinning* const _ccthin, /*optional*/
                         const vector<point>* const _skel_vts,
                         const vector<ivec2>* const _skel_edges,
                         const vector<uTriFace>* const _skel_faces,
                         vector<float>& _site_scalar);
} // namespace apps
} // namespace voxelvoro