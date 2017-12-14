#include "highlevelalgo.h"

#include <chrono>

#include <TriMesh.h>
#include <TriMesh_algo.h>
#include <plyall.h>

#include "octree.h"
#include "voxelapps.h"
#include "..\include\highlevelalgo.h"

using trimesh::TriMesh;

namespace voxelvoro
{
	bool preprocessVoro( VoroInfo & _voro, const shared_ptr<Volume3DScalar>& _vol, bool _need_euler )
	{
		timer t_tagging, t_merging;

		std::cout << "Starting preprocessing voro..." << std::endl;
		t_tagging.start();
		if ( !_voro.tagVtsUsingUniformVol( _vol ) )
		{
			cout << "Error: couldn't tag voro vts." << endl;
			return false;
		}
		t_tagging.stop();

		eulerchar euler_struct;
		if ( _need_euler )
		{
			_voro.computeEulerChar( euler_struct );
			euler_struct.logToConsole( "pruned voro diagram" );
		}

		auto max_len = std::max( std::max( _vol->getSizeX(), _vol->getSizeY() ), _vol->getSizeZ() );
		float eps = _voro.getInsidePartSize() * ( 1.0 / max_len ) * 0.000001;

		t_merging.start();
		cout << "before merging voro size: "
			<< _voro.geom().numVts() << "/"
			<< _voro.geom().numEdges() << "/"
			<< _voro.geom().numFaces() << endl;
		_voro.mergeCloseVts( eps, true );
		cout << "after merging voro size: "
			<< _voro.geom().numVts() << "/"
			<< _voro.geom().numEdges() << "/"
			<< _voro.geom().numFaces() << endl;
		t_merging.stop();
		
		cout << "time -> VD in/out tagging: " << t_tagging.elapseMilli().count() << " ms" << endl;
		cout << "time -> IVD-merge: " << t_merging.elapseMilli().count() << " ms" << endl;
		if ( _need_euler )
		{	// check euler numbers after merging identical vts
			_voro.computeEulerChar( euler_struct );
			euler_struct.logToConsole(
				( string( "after merging vts within " ) + std::to_string( eps ) ).c_str() );
		}
		std::cout << "Done preprocessing voro." << std::endl;
		return true;
	}
	
	ExportErrCode exportScalarFieldOnVoxelSurface(
		const char* _mc_geom_filename,
		const char* _mc_msure_filename,
		const char* _mc_order_filename,
		const char* _skel_filename, // optional
		const char* _output_filename,
		const char* _mc_msure_name,
		float _smooth_r,
		VoroInfo& _voro )
	{
		vector<point> mc_vts;
		vector<float> mc_msure;
		vector<int> mc_order;
		vector<point> skel_vts;
		vector<ivec2> skel_edges;
		vector<uTriFace> skel_faces;
		if ( readMedialCurveInfo( _mc_geom_filename, _mc_msure_filename, _mc_order_filename, _skel_filename,
			mc_vts, mc_msure, mc_order, &skel_vts, &skel_edges, &skel_faces ) != ImportErrCode::SUCCESS )
		{
			std::cout << "Failed reading medial-curve information!" << std::endl;
			return ExportErrCode::FAILURE;
		}
		std::cout << "Done: reading medial-curve information." << std::endl;

		vector<float> scalar_on_voxel_vts;
		auto msure_name = string( _mc_msure_name );
		voxelvoro::apps::MCMeasure msure_type;
		if ( msure_name == "bt3" )
		{
			msure_type = voxelvoro::apps::MCMeasure::BT3;
		}
		else if ( msure_name == "bt2" )
		{
			msure_type = voxelvoro::apps::MCMeasure::BT3;
		}
		else if ( msure_name == "traveldist" )
		{
			msure_type = voxelvoro::apps::MCMeasure::RIDGE;
		}
		else if ( msure_name == "length" )
		{
			msure_type = voxelvoro::apps::MCMeasure::LENGTH;
		}
		else if ( msure_name == "seglabel" )
		{
			msure_type = voxelvoro::apps::MCMeasure::SegmentLabel;
		}
		voxelvoro::apps::assignScalarToSites( _voro, 
			mc_vts, mc_msure, mc_order, msure_type, _smooth_r, 
			nullptr, 
			skel_vts.empty() ? nullptr : &skel_vts, 
			skel_vts.empty() ? nullptr : &skel_edges, 
			skel_vts.empty() ? nullptr : &skel_faces,
			scalar_on_voxel_vts );
		std::cout << "Done: assigning scalar field on voxel vts." << std::endl;

		//writeScalars( _output_filename, scalar_on_voxel_vts );
		std::ofstream scalar_os( _output_filename );
		if ( !scalar_os.is_open() )
		{
			std::cout << "Couldn't write to file " << _output_filename << std::endl;
			return ExportErrCode::OUTPUT_NOT_OPEN;
		}
		scalar_os << "# Each scalar corresponds to a voxel vertex. "
			<< "This file should pair with the.node file containing voxel vts. " << std::endl;
		scalar_os << scalar_on_voxel_vts.size() << endl;
		for ( auto s : scalar_on_voxel_vts )
			scalar_os << s << std::endl;
		scalar_os.close();
		std::cout << "Done: outputting scalar field to file " << _output_filename << std::endl;

		return ExportErrCode::SUCCESS;
	}
	
	ExportErrCode exportInsideVoroMesh(
		VoroInfo& _voro,
		const shared_ptr<Volume3DScalar>& _vol,
		const char* _mesh_file,
		bool _need_euler,
		bool _collapse_degenerate_edges,
		bool _inside_only, bool _finite_only,
		double _tt/*thinning threshold*/ )
	{
		timer t_thin, t_write_IVD, t_write_thinIVD;

		if ( !preprocessVoro(_voro, _vol, _need_euler) )
		{
			cout << "Error processing voro!" << endl;
			return ExportErrCode::FAILURE;
		}

		/*if ( _collapse_degenerate_edges )
		{
		shared_ptr<vector<point>> vts = make_shared<vector<point>>( mesh.vertices );
		shared_ptr<vector<uTriFace>> faces = make_shared<vector<uTriFace>>();
		std::for_each( mesh.faces.begin(), mesh.faces.end(), [ & ]( const TriMesh::Face& _f ) {
		faces->push_back( uTriFace(_f.v) ); }
		);
		shared_ptr<voxelvoro::EdgeCollapser> ec = make_shared<voxelvoro::TopoPreservEdgeCollapser>( vts, faces );
		cout << "edge collapser pre-processing ..." << endl;
		ec->preProcess();
		cout << "edge collapser processing done!" << endl;
		cout << "edge collapser collapsing degenerate edges ..." << endl;
		ec->collapse();
		cout << "edge collapser collapsing done!" << endl;

		eulerchar euler_struct;
		util::computeEulerChar( *faces, euler_struct );
		euler_struct.logToConsole( "voro after degenerate edges collapsed" );

		output_tri_faces.clear();
		output_tri_faces.insert( output_tri_faces.end(), faces->begin(), faces->end() );
		mesh.faces = output_tri_faces;
		trimesh::remove_unused_vertices( &mesh );
		}*/

		if ( std::string( _mesh_file ).find( ".off" ) != std::string::npos )
		{
			vector<uTriFace> output_tri_faces;
			vector<int> valid_face_indices;
			_inside_only ? _voro.getValidFaces( valid_face_indices ) :
				_finite_only ? _voro.getFiniteFaces( valid_face_indices ) : 
				_voro.getAllFaces( valid_face_indices );
			cout << "# faces about to write to file: " << valid_face_indices.size() << endl;
			_voro.getFaces( valid_face_indices, output_tri_faces );

			TriMesh mesh;
			mesh.vertices = _voro.geom().getAllVts();
			
			// convert vertices to original mesh space
			point v_trans;
			for ( auto& v : mesh.vertices )
			{
				SpaceConverter::fromVoxToModel( v, v_trans, _vol->getVoxToModelMat() );
				v = v_trans;
			}

			mesh.faces.clear();
			mesh.faces.reserve( output_tri_faces.size() );
			for ( auto fi = 0; fi < output_tri_faces.size(); ++fi )
			{
				const auto& f = output_tri_faces[ fi ];
				mesh.faces.push_back( ( TriMesh::Face )f.data() );
			}
			trimesh::remove_unused_vertices( &mesh );
			mesh.write( _mesh_file );
		}// else branch: output in .off format
		else if ( std::string( _mesh_file ).find( ".ply" ) != std::string::npos )
		{
			vector<point> vts;
			vector<ivec2> edges;
			vector<uTriFace> tri_faces;
			vector<float> vts_msure, edges_msure, faces_msure;
			_voro.extractInsideWithMeasure( MeasureForMA::LAMBDA,
				vts, edges, tri_faces,
				vts_msure, edges_msure, faces_msure );
			std::cout << "inside part extracted. " << std::endl;
			if ( vts.empty() )
			{
				cout << "IVD: 0 vts; IO skipped." << endl;
				goto END_PLY_BRANCH;
			}
			{ // local scope to allow goto to skip
				// convert vertices to original mesh space
				point v_trans;
				auto vts_trans = vts;
				for ( auto& v : vts_trans )
				{
					SpaceConverter::fromVoxToModel( v, v_trans, _vol->getVoxToModelMat() );
					v = v_trans;
				}

				// output info to file
				t_write_IVD.start();
				writeToPLY( _mesh_file, vts_trans, edges, tri_faces,
					vts_msure, edges_msure, faces_msure );
				t_write_IVD.stop();
				auto min_max_msure = std::minmax_element( faces_msure.begin(), faces_msure.end() );
				cout << "range of measure: " << "[" << *min_max_msure.first << "," << *min_max_msure.second << "]" << endl;

				float eps = _voro.getInsidePartSize() * 0.000001f;
				if ( !util::is_equal( _tt, 0.0, (double)eps ) && _tt > 0.0 )
				{
					// output pruned mesh
					t_thin.start();

					CellComplexThinning ccthin;
					cellcomplex cc( vts_trans, edges, tri_faces );
					ccthin.setup( &cc );
					ccthin.assignElementValues( vts_msure, edges_msure, faces_msure );
					ccthin.preprocess();
					auto minmax_fmsure = std::minmax_element( faces_msure.begin(), faces_msure.end() );
					//float thresh = ( *minmax_fmsure.second ) * 0.1f;
					float thresh = _voro.getInsidePartSize() * _tt;
					ccthin.prune( thresh, thresh, false );
					ccthin.remainingCC( vts_trans, edges, tri_faces );

					t_thin.stop();

					if ( _need_euler )
					{
						eulerchar euler_struct;
						util::computeEulerChar( vts_trans, edges, tri_faces, euler_struct );
						euler_struct.logToConsole( "voro diagram after thinning" );
					}
					if ( vts_trans.empty() )
					{
						cout << "Thinned IVD: 0 vts remaining; IO skipped." << endl;
						goto END_PLY_BRANCH;
					}

					// output the thinned result to file
					vts_msure.clear(); edges_msure.clear(); faces_msure.clear();
					std::string _thin_file = _mesh_file;
					auto end_str = std::string( "_thinned" ) + std::to_string( _tt );
					_thin_file.insert( _thin_file.find( ".ply" ), end_str );
					cout << "writing remaining cc to file: " << _thin_file << endl;
					t_write_thinIVD.start();
					writeToPLY( _thin_file.c_str(), vts_trans, edges, tri_faces,
						vts_msure, edges_msure, faces_msure );
					t_write_thinIVD.stop();
				}
			} // local scope for goto to skip
		END_PLY_BRANCH:
			cout << "time -> thinning: " << t_thin.elapseMilli().count() << " ms" << endl;
			cout << "time -> write IVD (I/O): " << t_write_IVD.elapseMilli().count() << " ms" << endl;
			cout << "time -> write thinned IVD (I/O): " << t_write_thinIVD.elapseMilli().count() << " ms" << endl;
		} // else branch: output in .ply format
		else
			return ExportErrCode::OUTPUT_NOT_SUPPORTED;

		return ExportErrCode::SUCCESS;
	}
	
	bool denseToSparse( const char * _src, const char * _dest )
	{
		shared_ptr<Volume3DScalar> src_vol;
		if ( voxelvoro::readMRC( _src, src_vol ) != voxelvoro::ImportErrCode::SUCCESS )
		{
			cout << "Error: cannot read volume from file: " << _src << endl;
			return false;
		}
		cout << "Constructing octree from dense volume..." << endl;
		OctreeVolume oct_vol( src_vol.get() );
		cout << "Done octree construction." << endl;
		cout << "Writing octree to file..." << endl;
		oct_vol.writeToFile( _dest );
		cout << "Done octree writting." << endl;
		return true;
	}
	int nClosedComponents( const cellcomplex & _cc )
	{
		int n_cc = 0;
		CellComplexThinning ccthin;
		//TODO: maybe change the argument type to const *?
		ccthin.setup( const_cast<cellcomplex*>( &_cc ) ); 
		vector<float> v_msure, e_msure, f_msure;
		// we want to remove everything except closed pockets. so actual measure values are irrelevant. 
		// we set measures all to 0's, and threshold to 1
		for ( auto i = 0; i < _cc.numVts(); ++i )
			v_msure.push_back( 0.0f );
		for ( auto i = 0; i < _cc.numEdges(); ++i )
			e_msure.push_back( 0.0f );
		for ( auto i = 0; i < _cc.numFaces(); ++i )
			f_msure.push_back( 0.0f );
		ccthin.assignElementValues( v_msure, e_msure, f_msure );
		// process and prune
		ccthin.preprocess();
		ccthin.prune( 1.0f, 1.0f, false );
		// get remaining cc
		auto remain_cc = ccthin.remainingCC();
		// num of closed components?
		// only worth counting when there are faces left
		if ( remain_cc.numFaces() == 0 )
			return 0;
		n_cc = remain_cc.compute_conn_cmpnts();
		return n_cc;
	}
}