#ifndef EXPORTER_H
#define EXPORTER_H

#include <vector>
using std::vector;

#include <reader.h> // tao's MRC reader

#include "Volume3DScalar.h"
#include "commondefs.h"
#include "voroinfo.h"

namespace voxelvoro
{
enum class ExportErrCode
{
    FAILURE,
    SUCCESS,
    INPUT_NOT_OPEN,      // cannot open input file
    OUTPUT_NOT_OPEN,     // cannot open output file
    INPUT_NOT_SUPPORTED, // input file extension not supported
    OUTPUT_NOT_SUPPORTED // output file extension not supported
};

// write the volume in a MRC file to a mathematica friendly file
int writeToMathematicaFromMRC(const shared_ptr<Volume3DScalar>& _vol,
                              const char* _outfile_name);
//
// write the boundary elements (vts & faces) into a mesh .off file from the
// given volume
// @note: the volume could be a dense or sparse rep.
ExportErrCode writeVolumeAsBoundaryMesh(const shared_ptr<Volume3DScalar>& _vol,
                                        const char* _outfile_base,
                                        bool _need_euler,
                                        bool _write_node = false);
//
// write the boundary vertices into a mesh .off file & a tetgen .node file from
// the given volume
ExportErrCode writeVolumeAsBoundaryPts(const shared_ptr<Volume3DScalar>& _vol,
                                       const char* _outfile_base,
                                       bool _write_node = false);
//
// write the given list of vertices to an output file (tetgen's .node file)
ExportErrCode writeToTetnodes(const vector<point>& _vts,
                              const char* _outfile_name);
//
// write given geometry & radius-per-vertex to the given .ma file (qmat file
// format)
ExportErrCode writeToDotma(const string& _dotma_filename,
                           const vector<point>& _vts,
                           const vector<ivec2> _all_edges,
                           const vector<uTriFace>& _tri_faces,
                           const vector<float>& _radii);
//
// write the inside part of a voro-diagram to .ma file (input format of QMAT)
ExportErrCode writeInsideVoroToDotMA(const VoroInfo& _voro,
                                     const char* _dotma_file);
//
// write vts, edges, and faces to a ply mesh file optionally with a measure for
// each element. user can also specify whether to output given site information
// assoc.ed with each face to the file
ExportErrCode
writeToPLY(const char* _ply_filename, const vector<point>& _output_vts,
           const vector<ivec2>& _output_edges,
           const vector<uTriFace>& _output_tris,
           const vector<float>& _vts_msure, const vector<float>& _edges_msure,
           const vector<float>& _faces_msure, bool _write_sites = false,
           const vector<point>* _sites = nullptr,
           const vector<ivec2>* _face_sites_ids = nullptr);
//
// write the radii filed out to a file (.r format)
// optionally transform radii using the given matrix
ExportErrCode writeRadiiToFile(const VoroInfo& _voro, const char* _r_file,
                               const trimesh::xform& _mat,
                               const vector<int>* _ids_ptr = nullptr);
//
// read in the medial axis (.off file) and the boundary vertices
// of the corresponding shape, estimate the radii field over the medial axis
int estimateRadiiField(const char* _ma_file_name,
                       const char* _bndry_pts_file_name,
                       const char* _radii_file_name);
//
// write a list of scalars to file
// ExportErrCode writeScalars( const char* _output_filename, vector<float>&
// _scalar_on_voxel_vts );
} // namespace voxelvoro

#endif
