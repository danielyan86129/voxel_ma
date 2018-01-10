#ifndef IMPORTERS_H
#define IMPORTERS_H

#include "commondefs.h"
#include "voroinfo.h"
#include "Volume3DScalar.h"

namespace voxelvoro
{
	enum class ImportErrCode
	{
		FAILURE,
		SUCCESS,
		// cannot open input file
		NOOPENINPUT,
		// doesn't support input volume file
		INVALID_VOL_FILE,
		// general unknown format error
		INVALID_FORMAT
	};

	//
	// read in a volume file (user should use this as much as possible instead of the other specific readers)
	ImportErrCode readVolume( const char* _vol_file, shared_ptr<Volume3DScalar>& _vol );
	//
	// read in an mrc volume
	ImportErrCode readMRC( const char* _mrc_file, shared_ptr<Volume3DScalar>& _vol );
	//
	// read an octree volume from a file (currently support .sof, .sog)
	ImportErrCode readOctree( const char* _filename, shared_ptr<Volume3DScalar>& _oct_vol );
	//
	// read in a list of points from .node file ( format used by tetgen )
	ImportErrCode readPts( const char* _node_file, vector<point>& _pts );
	//
	// import voro info from a set of tetgen files (deduced from given file base name)
	// @param _sites_pos_file: optionally import the contributing sites from the given file
	ImportErrCode readVoroInfo( const char* _basename, VoroInfo& _voro, 
		bool _need_euler, const char* _sites_pos_file = nullptr, shared_ptr<Volume3DScalar> _vol = nullptr );
	//
	// import mesh from the given file name (file formats: .ply | .off)
	// the file could describe a manifold/non-manifold 2-cell complex
	// returns a cellcomplex representing the mesh
	ImportErrCode readMesh( const string& _filename, cellcomplex& _cc, bool _finalize_cc = true );
	//
	// read vts, edges, and faces from a ply mesh file along with measures
	ImportErrCode readFromPLY( const char* _ply_filename,
		vector<point>& _output_vts, vector<ivec2>& _output_edges, vector<uTriFace>& _output_tris,
		vector<float>& _vts_msure, vector<float>& _edges_msure, vector<float>& _faces_msure );
	//
	// read in a list of floats from a file
	ImportErrCode readNumberList( const char* _nums_filename, vector<float>& _nums );
	//
	// read in a few pieces of information describing a medial-curve structure
	// from a bunch of files
	ImportErrCode readMedialCurveInfo(
		const char * _mc_geom_filename,
		const char * _mc_msure_filename,
		const char * _mc_order_filename,
		const char* _skel_name,
		vector<point>& _mc_vts, vector<float>& _mc_msure, vector<int>& _mc_order, 
		vector<point>* _skel_vts, vector<ivec2>* _skel_edges, vector<uTriFace>* _skel_faces );
}

#endif