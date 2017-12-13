#ifndef EXPORTER_H
#define EXPORTER_H

#include <vector>
using std::vector;

#include <reader.h> // tao's MRC reader

#include "commondefs.h"
#include "voroinfo.h"
#include "Volume3DScalar.h"

namespace voxelvoro
{
	enum class ExportErrCode
	{
		FAILURE,
		SUCCESS,
		INPUT_NOT_OPEN, // cannot open input file
		OUTPUT_NOT_OPEN, // cannot open output file
		INPUT_NOT_SUPPORTED, // input file extension not supported
		OUTPUT_NOT_SUPPORTED // output file extension not supported
	};

	// write the volume in a MRC file to a mathematica friendly file
	int writeToMathematicaFromMRC( const shared_ptr<Volume3DScalar>& _vol, const char * _outfile_name );
	//
	// write the boundary elements (vts & faces) into a mesh .off file from the given volume
	// @note: the volume could be a dense or sparse rep. 
	ExportErrCode writeVolumeAsBoundaryMesh( const shared_ptr<Volume3DScalar>& _vol, const char * _outfile_base, bool _need_euler );
	//
	// write the boundary vertices into a mesh .off file & a tetgen .node file from the given volume
	ExportErrCode writeVolumeAsBoundaryPts( const shared_ptr<Volume3DScalar>& _vol, const char * _outfile_base );
	//
	// write the given list of vertices to an output file (tetgen's .node file)
	ExportErrCode writeToTetnodes( const vector<point>& _vts, const char* _outfile_name );
	//
	// write the inside part of a voro-diagram to .ma file (input format of QMAT)
	ExportErrCode writeInsideVoroToDotMA( const VoroInfo& _voro, const char* _dotma_file );
	//
	// write vts, edges, and faces to a ply mesh file along with measures
	ExportErrCode writeToPLY( const char* _ply_filename,
		const vector<point>& _output_vts, const vector<ivec2>& _output_edges, const vector<uTriFace>& _output_tris,
		const vector<float>& _vts_msure, const vector<float>& _edges_msure, const vector<float>& _faces_msure );
	//
	// write the radii filed out to a file (.r format)
	// optionally transform radii using the given matrix
	ExportErrCode writeRadiiToFile( const VoroInfo& _voro, const char* _r_file, const trimesh::xform& _mat );
	//
	// read in the medial axis (.off file) and the boundary vertices 
	// of the corresponding shape, estimate the radii field over the medial axis
	int estimateRadiiField( const char* _ma_file_name, const char * _bndry_pts_file_name, 
		const char* _radii_file_name);
	//
	// write a list of scalars to file
	//ExportErrCode writeScalars( const char* _output_filename, vector<float>& _scalar_on_voxel_vts );
}

#endif
