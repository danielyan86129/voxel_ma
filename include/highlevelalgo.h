#pragma once
#include "importers.h"
#include "exporters.h"

namespace voxelvoro
{
	//
	// process voro diagram to properly tag its element (in-out)
	// and simply degenerate elements
	bool preprocessVoro( 
		VoroInfo & _voro, 
		const shared_ptr<Volume3DScalar>& _vol, 
		bool _need_euler );
	//
	// Read medial-curve related information from a set of specified input files,
	// assign a scalar field using that information,
	// and finally output the scalar field to the specified output file.
	ExportErrCode exportScalarFieldOnVoxelSurface(
		const char* _mc_geom_filename,
		const char* _mc_msure_filename,
		const char* _mc_order_filename,
		const char* _skel_filename, // optional
		const char* _output_filename,
		const char* _mc_msure_type,
		float _smooth_r,
		VoroInfo& _voro );
	//
	// extract and export the part of the voro-info object that's "inside" the volume
	ExportErrCode exportInsideVoroMesh( 
		VoroInfo& _voro, 
		const shared_ptr<Volume3DScalar>& _vol, 
		const char* _mesh_file,
		bool _need_euler, 
		bool _collapse_degenerate_edges, 
		bool _inside_only, bool _finite_only,
		double _tt/*thinning threshold*/ );
	//
	// convert a dense volume (stored in .mrc file) to a sparse volume (stored in .sof file)
	bool denseToSparse( const char* _src, const char* _dest );
	//
	// return the num of maximally closed sub-complex in the given cell complex
	int nClosedComponents( const cellcomplex& _cc );
}