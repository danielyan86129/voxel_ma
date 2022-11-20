#pragma once
#include "exporters.h"
#include "importers.h"
#include <tetgen/tetgen.h>

namespace voxelvoro
{
/**************************************
**  Tetgen structures related functions
****************************************/

/**
 * Commented out due to incompatible with tetgen1.6
 */
// //
// // my version of making a tetgenmesh from a tetgenio
// void my_tetrahedralize(char* switches, tetgenio* _in, tetgenio* out,
//                        tetgenmesh& m, tetgenio* addin = nullptr,
//                        tetgenio* bgmin = nullptr);
// void my_tetrahedralize(tetgenbehavior* _b, tetgenio* _in, tetgenio* _out,
//                        tetgenmesh& _m, tetgenio* _addin = nullptr,
//                        tetgenio* _bgmin = nullptr);
//
// convert a list of points to tetgen's rep
// NOTE:
//   No mem check. E.g. if @param _tetio.numofpoints is big enough,
//   we assume there is enough mem for the points
void pts2tetgen(const vector<point> _P, tetgenio& _tetio);
// compute VD of given voxel shape
void computeVD(const shared_ptr<Volume3DScalar>& _vol, VoroInfo& _voro);
//
// compute VD (or VD core) given sites:
//   tetrahedralize -> dual (voro) -> (optionally) extract core
void computeVD(tetgenio& _tet_in, VoroInfo& _voro,
               shared_ptr<Volume3DScalar> _vol = nullptr);
/*void computeVD( const vector<point>& _sites,
        VoroInfo& _voro,
        shared_ptr<Volume3DScalar> _vol = nullptr );*/
//
// process voro diagram to properly tag its element (in-out)
// and simply degenerate elements
bool preprocessVoro(VoroInfo& _voro, const shared_ptr<Volume3DScalar>& _vol,
                    bool _need_euler);
//
// Read medial-curve related information from a set of specified input files,
// assign a scalar field using that information,
// and finally output the scalar field to the specified output file.
ExportErrCode exportScalarFieldOnVoxelSurface(
    const char* _mc_geom_filename, const char* _mc_msure_filename,
    const char* _mc_order_filename,
    const char* _skel_filename, // optional
    const char* _output_filename, const char* _mc_msure_type, float _smooth_r,
    const VoroInfo& _voro);
//
// extract and export the part of the voro-info object that's "inside" the
// volume
ExportErrCode
exportInsideVoroMesh(VoroInfo& _voro, const shared_ptr<Volume3DScalar>& _vol,
                     const char* _mesh_file, bool _need_euler,
                     bool _collapse_degenerate_edges, bool _inside_only,
                     bool _finite_only,
                     int _full_or_thin, /*0: thin, 1: full, 2: both*/
                     bool _write_attrib, bool _out_qmat,
                     const vector<float>& _tt /*thinning threshold(s)*/);
//
// convert a dense volume (stored in .mrc file) to a sparse volume (stored in
// .sof file)
bool denseToSparse(const char* _src, const char* _dest);
//
// return the num of maximally closed sub-complex in the given cell complex
int nClosedComponents(const cellcomplex& _cc);
} // namespace voxelvoro