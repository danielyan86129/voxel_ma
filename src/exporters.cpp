#include "exporters.h"

#include <fstream>
#include <iostream>
#include <sstream> 
#include <memory>
#include <vector>
#include <string>
#include <cassert>
#include <chrono>

#include <TriMesh.h>
#include <TriMesh_algo.h>
#include <plyall.h>

#include "surfacing.h"
#include "geomalgo.h"
#include "edgecollapse.h"
#include "ccthin.h"
#include "importers.h"

using std::cout; using std::endl; using std::vector;
using std::shared_ptr; 
using std::string; 
using std::ifstream; using std::ofstream;
using namespace trimesh;

namespace voxelvoro
{
	int writeToMathematicaFromMRC( const shared_ptr<Volume3DScalar>& _vol, const char * _outfile_name )
	{
		int ret = 0;

		cout << "Start writing uniform volume to output file "<< _outfile_name << endl;
		std::stringstream outstr;
		outstr << "{ (*beginning of everything*)"<<endl;

		outstr << "{ (*begin: volume*)"<<endl;
		for (int i = 0; i < _vol->getSizeX(); ++i)
		{
			outstr << "{";
			
			for (int j = 0; j < _vol->getSizeX(); ++j)
			{
				outstr << "{";
				for (int k = 0; k < _vol->getSizeX(); ++k)
				{
					auto val = _vol->getDataAt(i,j,k);
					if ( val < 0 )
						val = 0;
					else
						val = 1;
					outstr << val << ",";
				}
				outstr.seekp(-1, outstr.cur);
				outstr << "},"<<endl;
			}
			
			// erase ", endl"
			outstr.seekp(-2, outstr.cur);
			outstr << endl;

			outstr << "},"<<endl;
		}
		// erase ", endl"
		outstr.seekp(-2, outstr.cur);
		outstr << endl;
		outstr << "} (*end: volume*)"<<endl;

		outstr << "} (*end of everything*)"<<endl;

		ofstream outfile(_outfile_name);
		if ( !outfile.is_open() )
		{
			cout << "Error: Cannot open file " << _outfile_name << endl;
			ret = -1;
			goto RETURN;
		}

		outfile << outstr.rdbuf();

		cout << "Done writing uniform volume to output file!" << endl;

	RETURN:
		outfile.close();
		return ret;
	}

	ExportErrCode writeVolumeAsBoundaryMesh( const shared_ptr<Volume3DScalar>& _vol, const char * _outfile_base, 
		bool _need_euler, bool _write_node )
	{
		timer t_extract_bndry, t_write_bndry;
		Surfacer surf;
		vector<point> bndry_vts;
		vector<uTriFace> bndry_faces;
		eulerchar euler;
		t_extract_bndry.start();
		auto retcode = surf.extractBoundaryVtsAndFaces( _vol, bndry_vts, bndry_faces, _need_euler ? &euler : nullptr );
		t_extract_bndry.stop();
		if ( retcode == SurfacerErrCode::INVALID_VOL_TYPE )
		{
			cout << "Error: volume type not supported!" << endl;
		} 
		else if ( retcode == SurfacerErrCode::EMPTY_VOL )
		{
			cout << "Error: empty volume was given!" << endl;
		}
		if ( retcode != SurfacerErrCode::SUCCESS )
		{
			return ExportErrCode::FAILURE;
		}
		if ( _need_euler )
		{
			euler.logToConsole( "boundary voxel mesh" );
		}

		auto mesh = trimesh::TriMesh();
		point v_trans;
		for ( const auto& v : bndry_vts )
		{
			v_trans = _vol->getVoxToModelMat() * v;
			mesh.vertices.push_back( v_trans );
		}
		for (const auto& f: bndry_faces )
		{
			mesh.faces.push_back( trimesh::TriMesh::Face( f[ 0 ], f[ 1 ], f[ 2 ] ) );
		}

		t_write_bndry.start();
		auto mesh_success = mesh.write( std::string( _outfile_base ) + ".off" );
		auto nodes_err = ExportErrCode::SUCCESS;
		if (_write_node )
			nodes_err = writeToTetnodes( bndry_vts, (std::string(_outfile_base) + ".node").c_str() );
		t_write_bndry.stop();
		cout << "time(I/O) -> extract boundary pts: " << t_extract_bndry.elapseMilli().count() << " ms" << endl;
		cout << "time(I/O) -> write boundary pts: " << t_write_bndry.elapseMilli().count() << " ms" << endl;
		if ( mesh_success && ( nodes_err == ExportErrCode::SUCCESS) )
		{
			return ExportErrCode::SUCCESS;
		}
		else
		{
			return mesh_success ? nodes_err : ExportErrCode::FAILURE;
		}
	}

	ExportErrCode writeVolumeAsBoundaryPts( const shared_ptr<Volume3DScalar>& _vol, const char * _outfile_base, bool _write_node )
	{
		timer t_extract_bndry, t_write_bndry;
		Surfacer surf;
		vector<point> bndry_vts;
		t_extract_bndry.start();
		auto retcode = surf.extractBoundaryVts( _vol, bndry_vts );
		t_extract_bndry.stop();
		if ( retcode == SurfacerErrCode::INVALID_VOL_TYPE )
		{
			cout << "Error: volume type not supported!" << endl;
		}
		else if ( retcode == SurfacerErrCode::EMPTY_VOL )
		{
			cout << "Error: given volume is empty!" << endl;
		}
		if ( retcode != SurfacerErrCode::SUCCESS )
		{
			return ExportErrCode::FAILURE;
		}

		// write a copy of boundary vts in the model space (useful for figure-making)
		trimesh::TriMesh mesh;
		point v_trans;
		for ( const auto& v : bndry_vts )
		{
			v_trans = _vol->getVoxToModelMat() * v;
			mesh.vertices.push_back( v_trans );
		}

		auto meshfile = std::string( _outfile_base ) + ".off";
		auto mesh_retcode = mesh.write( meshfile );
		if ( mesh_retcode == false )
		{
			cout << "failed to write boundary vts to \"" << meshfile << "\"..." << endl;
		}

		// output vts to .node file in voxel space (so that tetgen result is more robust)
		if (_write_node)
		{
			t_write_bndry.start();
			auto nodefile = ( std::string( _outfile_base ) + ".node" );
			cout << "writing boundary vts to \"" << nodefile << "\"..." << endl;
			auto nodes_retcode = writeToTetnodes( bndry_vts, nodefile.c_str() );
			if ( nodes_retcode != ExportErrCode::SUCCESS )
			{
				cout << "Error: failed to write boundary vts to .node file." << endl;
				return nodes_retcode;
			}
		}
		t_write_bndry.stop();

		cout << "time -> extract boundary pts: " << t_extract_bndry.elapseMilli().count() <<" ms"<< endl;
		cout << "time(I/O) -> write boundary pts: " << t_write_bndry.elapseMilli().count() << " ms" << endl;
		return mesh_retcode ? ExportErrCode::SUCCESS : ExportErrCode::FAILURE;
	}

	ExportErrCode writeToTetnodes( const vector<point>& _vts, const char * _outfile_name )
	{
		ofstream out_file( _outfile_name );
		if ( !out_file.is_open() )
		{
			cout << "Error: couldn't open " << _outfile_name << endl;
			return ExportErrCode::OUTPUT_NOT_OPEN;
		}

		out_file << _vts.size() << " " << 3 << " " << 0 << " " << 0 << endl;
		for ( int i = 0; i < _vts.size(); ++i )
		{
			const auto& v = _vts[ i ];
			out_file << i << " " << v[ 0 ] << " " << v[ 1 ] << " " << v[ 2 ] << endl;
		}

		out_file.close();

		return ExportErrCode::SUCCESS;
	}

	ExportErrCode writeToDotma( const string& _dotma_filename,
		const vector<point>& _vts, const vector<ivec2> _all_edges, const vector<uTriFace>& _tri_faces,
		const vector<float>& _radii )
	{
		ofstream ofile( _dotma_filename );
		if ( !ofile.is_open() )
		{
			cout << "Error: cannot open " << _dotma_filename << endl;
			return ExportErrCode::FAILURE;
		}

		// write out the file's header (some simple stats)
		ofile
			<< _vts.size() << " "
			<< _all_edges.size() << " "
			<< _tri_faces.size() << endl;

		// write out all vts & associated radius
		float r_min = std::numeric_limits<float>::max();
		float r_max = -999.0f;
		for ( auto i  = 0; i < _vts.size(); ++i )
		{
			const auto& v = _vts[i];
			auto r = _radii[ i ];
			r_min = std::min( r, r_min );
			r_max = std::max( r, r_max );
			ofile << "v " << v[ 0 ] << " " << v[ 1 ] << " " << v[ 2 ] << " " << r << endl;
		}
		//cout << "range of radius: [" << r_min << ", " << r_max << "]" << endl;

		// write out all edges
		for ( auto it = _all_edges.begin(); it != _all_edges.end(); ++it )
		{
			ofile << "e " << ( *it )[ 0 ] << " " << ( *it )[ 1 ] << endl;
		}

		// write out all faces
		for ( const auto& t : _tri_faces )
		{
			ofile << "f " << t[ 0 ] << " " << t[ 1 ] << " " << t[ 2 ] << endl;
		}

		ofile.close();

		return ExportErrCode::SUCCESS;
	}
	ExportErrCode writeInsideVoroToDotMA( const VoroInfo & _voro, const char * _dotma_filename )
	{
		vector<int> vts_indices;
		vector<ivec2> total_tri_edges;
		vector<uTriFace> total_tri_faces;
		_voro.getPartInside( vts_indices, total_tri_edges, total_tri_faces );
		const auto& voro_geom = _voro.geom();

		ofstream ofile( _dotma_filename );
		if ( !ofile.is_open() )
		{
			cout << "Error: cannot open " << _dotma_filename << endl;
			return ExportErrCode::FAILURE;
		}

		// write out the file's header (some simple stats)
		ofile 
			<< vts_indices.size() << " " 
			<< total_tri_edges.size() << " " 
			<< total_tri_faces.size() << endl;

		// write out all vts & associated radius
		auto radii = _voro.getRadii();
		float r_min = std::numeric_limits<float>::max();
		float r_max = -999.0f;
		for ( auto vi : vts_indices )
		{
			const auto& v = voro_geom.getVert( vi );
			auto r = radii[ vi ];
			r_min = std::min( r, r_min );
			r_max = std::max( r, r_max );
			ofile << "v " << v[ 0 ] << " " << v[ 1 ] << " " << v[ 2 ] << " " << r << endl;
		}
		//cout << "range of radius: [" << r_min << ", " << r_max << "]" << endl;

		// write out all edges
		for ( auto it = total_tri_edges.begin(); it != total_tri_edges.end(); ++it )
		{
			ofile << "e " << ( *it )[ 0 ] << " " << ( *it )[ 1 ] << endl;
		}

		// write out all faces
		for ( const auto& t : total_tri_faces )
		{
			ofile << "f " << t[ 0 ] << " " << t[ 1 ] << " " << t[ 2 ] << endl;
		}

		ofile.close();

		return ExportErrCode::SUCCESS;
	}

	ExportErrCode writeToPLY( const char* _ply_filename,
		const vector<point>& _output_vts, const vector<ivec2>& _output_edges, const vector<uTriFace>& _output_tris,
		const vector<float>& _vts_msure, const vector<float>& _edges_msure, const vector<float>& _faces_msure,
		bool _write_sites, const vector<point>* _sites, const vector<ivec2>* _face_sites_ids )
	{
		struct Vertex
		{
			float x; float y; float z;
			unsigned char r, g, b;
			float s; // measure
		};
		struct Edge
		{
			int v1; int v2;
			unsigned char r, g, b;
			float s; // measure
		};
		struct Face
		{
			unsigned char nvts;
			int verts[ 3 ];
			unsigned char sites_l; // len of sites array
			float sites[ 6 ]; // 2*3 coords for two sites positions
			float s; // measure
			unsigned char r, g, b; // pseudo-color of measure
		};
		std::map<std::string, PlyProperty> vert_props;
		vert_props[ "x" ] = { "x", Float32, Float32, offsetof( Vertex, x ), PLY_SCALAR, 0, 0, 0 };
		vert_props[ "y" ] = { "y", Float32, Float32, offsetof( Vertex, y ), PLY_SCALAR, 0, 0, 0 };
		vert_props[ "z" ] = { "z", Float32, Float32, offsetof( Vertex, z ), PLY_SCALAR, 0, 0, 0 };
		if ( !_vts_msure.empty() )
		{
			/*vert_props[ "red" ] = { "red", Uint8, Uint8, offsetof( Vertex, r ), PLY_SCALAR, 0, 0, 0 };
			vert_props[ "green" ] = { "green", Uint8, Uint8, offsetof( Vertex, g ), PLY_SCALAR, 0, 0, 0 };
			vert_props[ "blue" ] = { "blue", Uint8, Uint8, offsetof( Vertex, b ), PLY_SCALAR, 0, 0, 0 };*/
		}
		std::map<std::string, PlyProperty> edge_props;
		edge_props[ "vertex1" ] = { "vertex1", Int32, Int32, offsetof( Edge, v1 ), PLY_SCALAR, 0, 0, 0 };
		edge_props[ "vertex2" ] = { "vertex2", Int32, Int32, offsetof( Edge, v2 ), PLY_SCALAR, 0, 0, 0 };
		if ( !_edges_msure.empty() )
		{
			/*edge_props[ "red" ] = { "red", Uint8, Uint8, offsetof( Edge, r ), PLY_SCALAR, 0, 0, 0 };
			edge_props[ "green" ] = { "green", Uint8, Uint8, offsetof( Edge, g ), PLY_SCALAR, 0, 0, 0 };
			edge_props[ "blue" ] = { "blue", Uint8, Uint8, offsetof( Edge, b ), PLY_SCALAR, 0, 0, 0 };*/
		}
		std::map<std::string, PlyProperty> face_props;
		face_props[ "vertex_indices" ] = {
			"vertex_indices", Int32, Int32, offsetof( Face, verts ),
			PLY_LIST, Uint8, Uint8, offsetof( Face,nvts ) };
		if ( !_faces_msure.empty() )
		{
			face_props[ "msure" ] = { "msure", Float32, Float32, offsetof( Face, s ), PLY_SCALAR, 0,0,0 };
			/*face_props[ "red" ] = { "red", Uint8, Uint8, offsetof( Face, r ), PLY_SCALAR, 0, 0, 0 };
			face_props[ "green" ] = { "green", Uint8, Uint8, offsetof( Face, g ), PLY_SCALAR, 0, 0, 0 };
			face_props[ "blue" ] = { "blue", Uint8, Uint8, offsetof( Face, b ), PLY_SCALAR, 0, 0, 0 };*/
		}
		if ( _write_sites )
		{
			face_props[ "sites" ] = { "sites", Float32, Float32, offsetof( Face, sites ), PLY_LIST, Uint8, Uint8, offsetof( Face, sites_l ) };
		}
		vector<Vertex> output_vts( _output_vts.size() );
		vector<Edge> output_edges( _output_edges.size() );
		vector<Face> output_faces( _output_tris.size() );
		// fill in the data to output lists (face might contain color)
		auto min_max_msure = std::minmax_element( _vts_msure.begin(), _vts_msure.end() );
		for ( auto i = 0; i < output_vts.size(); ++i )
		{
			auto& o_p = output_vts[ i ];
			o_p.x = _output_vts[ i ][ 0 ];
			o_p.y = _output_vts[ i ][ 1 ];
			o_p.z = _output_vts[ i ][ 2 ];
			/*if ( !_vts_msure.empty() )
			{
			RGBColor c = (RGBColor)util::GetColour(
			_vts_msure[ i ], *min_max_msure.first, *min_max_msure.second );
			o_p.r = c.r();
			o_p.g = c.g();
			o_p.b = c.b();
			}*/
		}
		min_max_msure = std::minmax_element( _edges_msure.begin(), _edges_msure.end() );
		for ( auto i = 0; i < output_edges.size(); ++i )
		{
			auto& o_e = output_edges[ i ];
			o_e.v1 = _output_edges[ i ][ 0 ];
			o_e.v2 = _output_edges[ i ][ 1 ];
			/*if (!_edges_msure.empty() )
			{
			RGBColor c = (RGBColor)util::GetColour(
			_edges_msure[ i ], *min_max_msure.first, *min_max_msure.second );
			o_e.r = c.r();
			o_e.g = c.g();
			o_e.b = c.b();
			}*/
		}
		min_max_msure = std::minmax_element( _faces_msure.begin(), _faces_msure.end() );
		for ( auto i = 0; i < output_faces.size(); ++i )
		{
			auto& o_f = output_faces[ i ];
			o_f.nvts = 3;
			o_f.verts[ 0 ] = _output_tris[ i ][ 0 ];
			o_f.verts[ 1 ] = _output_tris[ i ][ 1 ];
			o_f.verts[ 2 ] = _output_tris[ i ][ 2 ];
			if ( !_faces_msure.empty() )
			{
				o_f.s = _faces_msure[ i ];
				/*RGBColor c = (RGBColor)util::GetColour(
					o_f.s, *min_max_msure.first, *min_max_msure.second );
				o_f.r = c.r();
				o_f.g = c.g();
				o_f.b = c.b();*/
			}
			if ( _write_sites )
			{
				o_f.sites_l = 6;
				auto site_ids = (*_face_sites_ids)[ i ];
				auto s1 = (*_sites)[ site_ids[ 0 ] ];
				auto s2 = ( *_sites )[ site_ids[ 1 ] ];
				memcpy( o_f.sites, s1.data(), sizeof( float ) * 3 );
				memcpy( &o_f.sites[3], s2.data(), sizeof( float ) * 3 );
			}
		}

		std::cout << "writing to ply file... " << std::endl;
		std::cout << "# vts/edges/faces: "
			<< output_vts.size() << "/"
			<< output_edges.size() << "/"
			<< output_faces.size() << std::endl;

		ply::PLYwriter ply_writer;
		if (
			ply_writer.write( _ply_filename, true, true, true,
				vert_props, edge_props, face_props,
				output_vts, output_edges, output_faces ) != ply::SUCCESS
			)
			return ExportErrCode::FAILURE;
		return ExportErrCode::SUCCESS;
	}

	ExportErrCode writeRadiiToFile( const VoroInfo& _voro, const char * _r_file, const trimesh::xform& _mat, const vector<int>* _ids_ptr )
	{
		if ( !_voro.radiiValid() )
		{
			cout << "Nothing written to file. Radii of the given voro is not correctly set yet!" << endl;
			return ExportErrCode::FAILURE;
		}

		ofstream rfile( _r_file );
		if ( !rfile.is_open() )
		{
			cout << "Error: cannot open file " << _r_file << endl;
			return ExportErrCode::FAILURE;
		}

		timer t_IO;
		t_IO.start();
		int cnt;
		if ( !_ids_ptr )
		{
			// write only valid vertices' radii
			cnt = 0;
			for ( size_t i = 0; i < _voro.geom().numVts(); ++i )
				if ( _voro.isVertexValid( i ) )
					cnt++;
			rfile << cnt << endl;

			const auto& radii = _voro.getRadii();
			auto orig_trans = _mat *point( 0, 0, 0 );
			for ( size_t i = 0; i < _voro.geom().numVts(); ++i )
			{
				if ( _voro.isVertexValid( i ) )
				{
					auto r = radii[ i ];
					// transform r
					r = trimesh::len( ( _mat * point( r, 0, 0 ) ) - orig_trans );
					auto v = _voro.geom().getVert( i );
					v = _mat * v;
					rfile << v[ 0 ] << " " << v[ 1 ] << " " << v[ 2 ] << " " << r << endl;
					cnt++;
				}
			}
		}
		else
		{
			// write given set of vertices
			const auto ids = *_ids_ptr;
			cnt = ids.size();
			rfile << cnt << endl;

			const auto& radii = _voro.getRadii();
			auto orig_trans = _mat *point( 0, 0, 0 );
			for ( size_t i : ids )
			{
				auto r = radii[ i ];
				// transform r
				r = trimesh::len( ( _mat * point( r, 0, 0 ) ) - orig_trans );
				auto v = _voro.geom().getVert( i );
				v = _mat * v;
				rfile << v[ 0 ] << " " << v[ 1 ] << " " << v[ 2 ] << " " << r << endl;
				cnt++;
			}
		}
		rfile.close();
		t_IO.stop();
		cout << "Done: " << cnt << " radii written to file " << _r_file << endl;
		cout << "time(I/O) -> write radii: " << t_IO.elapseMilli().count() << " ms" << endl;
		return ExportErrCode::SUCCESS;
	}

	int estimateRadiiField( const char* _ma_file_name, const char * _bndry_pts_file_name, 
		const char* _radii_file_name)
	{
		int ret = 0;
		float* bndry_pts, *radii;
		std::shared_ptr<trimesh::KDtree> tree;
		std::string ln;
		std::stringstream ss;
		ifstream bndry_file;
		ofstream radii_file;

		// read medial axis
		auto m_ptr = TriMesh::read(_ma_file_name);
		auto ma_sptr = shared_ptr<TriMesh>(m_ptr);
		if (!ma_sptr)
		{
			cout << "Error reading file " << _ma_file_name << endl;
			ret = -1;
			goto RETURN;
		}

		// read boundary vertices file
		//typedef trimesh::Vec<3, int> intpoint;
		bndry_file.open( _bndry_pts_file_name );
		if ( !bndry_file.is_open() )
		{
			cout << "Error: couldn't open file "<<_bndry_pts_file_name<<endl;
			ret = -1;
			goto RETURN;
		}

		// read file header
		int n_bndry_vts = 0;
		while ( !bndry_file.eof() )
		{
			ln.clear();
			std::getline(bndry_file, ln);
			// skip comment
			if ( ln[0] != '#' )
				break;
		}
		ss << ln;
		ss >> n_bndry_vts;

		bndry_pts = new float[n_bndry_vts * 3];

		int i = 0;
		int v_id;
		float x, y, z;
		while ( !bndry_file.eof() )
		{
			ln.clear();
			std::getline(bndry_file, ln);
			// skip comment
			if ( ln.empty() || ln[0] == '#' )
				continue;

			ss.str(ln);
			ss.clear();
			ss >> v_id >> x >> y >> z;
			bndry_pts[3*i+0] = x;
			bndry_pts[3*i+1] = y;
			bndry_pts[3*i+2] = z;
			i ++;
		}
		bndry_file.close();
		assert(i == n_bndry_vts);
		
		// now use kdtree to find nearest boundary point for each ma point
		// their distance is the radius for that ma point
		tree = shared_ptr<trimesh::KDtree>( 
			new trimesh::KDtree (bndry_pts, n_bndry_vts) 
			);
		radii = new float[ma_sptr->vertices.size()];
		ma_sptr->need_bbox();
		for (int i = 0; i < ma_sptr->vertices.size(); ++i)
		{
			auto v = ma_sptr->vertices[i];
			auto closest = tree->closest_to_pt( v.data(), ma_sptr->bbox.radius()*ma_sptr->bbox.radius() );
			
			radii[i] = trimesh::dist( trimesh::point(closest), v );
		}

		// write radii to output file
		radii_file.open(_radii_file_name);
		radii_file << ma_sptr->vertices.size() << endl;
		for ( int i = 0; i < ma_sptr->vertices.size(); ++i )
		{
			auto r = radii[i];
			auto v = ma_sptr->vertices[i];
			radii_file << v[0] << " " << v[1] << " " << v[2] << " " << r << endl;
		}
		radii_file.close();

		// clean-up
		delete[] bndry_pts;
		delete[] radii;

RETURN:
		return ret;
	}
}