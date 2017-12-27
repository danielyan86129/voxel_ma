#include "importers.h"
#include "densevolume.h"
#include "octree.h"
#include <reader.h> // Tao's volume reader
#include <sstream> 
#include <cassert>
#include <plyall.h>
#include <filesystem>

namespace voxelvoro {
	using std::cout;
	using std::endl;
	namespace fs = std::experimental::filesystem;

	ImportErrCode readVolume( const char * _vol_file, shared_ptr<Volume3DScalar>& _vol )
	{
		timer t_create_vol;
		t_create_vol.start();
		std::string filename( _vol_file );
		auto found = filename.find_last_of( '.' );
		auto ext = filename.substr( found + 1 );
		std::transform( ext.begin(), ext.end(), ext.begin(), ::tolower );
		if ( ext == "mrc" )
		{
			return readMRC( _vol_file, _vol );
		}
		else if ( ext == "sof" || ext == "sog" )
		{
			return readOctree( _vol_file, _vol );
		}
		else
		{
			return ImportErrCode::INVALID_VOL_FILE;
		}
		t_create_vol.stop();
		cout << "time -> create volume: " << t_create_vol.elapseMilli().count() << " ms" << endl;
	}

	ImportErrCode readMRC( const char * _mrc_file, shared_ptr<Volume3DScalar>& _dense_vol )
	{
		char infile_name[ 1024 ];
		strcpy( infile_name, _mrc_file );

		cout << "Start reading MRC file " << _mrc_file << endl;

		auto mrc_reader = shared_ptr<VolumeReader>( MRCReaderPicker::pick( infile_name ) );
		auto vol = shared_ptr<Volume>( mrc_reader->getVolume() );
		_dense_vol = shared_ptr<DenseVolume>( new DenseVolume( vol ) );
		cout << "Done reading MRC file!" << endl;

		return ImportErrCode::SUCCESS;
	}

	ImportErrCode readOctree( const char * _filename, shared_ptr<Volume3DScalar>& _oct_vol )
	{
		char infile_name[ 1024 ];
		strcpy( infile_name, _filename );

		cout << "Start reading .SOF file " << _filename << endl;

		// TODO: handle failure of file reading (issued from within octree construction)
		_oct_vol = std::static_pointer_cast<Volume3DScalar>(
			shared_ptr<OctreeVolume>( new OctreeVolume( infile_name ) )
			);
		if ( !_oct_vol )
			cout << "octree-to-parent-volume-type conversion failed!" << endl;
		return ImportErrCode::SUCCESS;
	}

	ImportErrCode readPts( const char * _node_file, vector<point>& _pts )
	{
		timer t_io;
		t_io.start();
		std::string ln;
		std::stringstream ss;
		ifstream bndry_file;

		// read boundary vertices file
		//typedef trimesh::Vec<3, int> intpoint;
		bndry_file.open( _node_file );
		if ( !bndry_file.is_open() )
		{
			cout << "Error: couldn't open file " << _node_file << endl;
			return ImportErrCode::NOOPENINPUT;
		}

		// read file header
		int n_bndry_vts = 0;
		while ( !bndry_file.eof() )
		{
			ln.clear();
			std::getline( bndry_file, ln );
			// skip comment
			if ( ln[ 0 ] != '#' )
				break;
		}
		ss << ln;
		ss >> n_bndry_vts;
		_pts.resize( n_bndry_vts );

		int i = 0;
		int v_id;
		float x, y, z;
		while ( !bndry_file.eof() )
		{
			ln.clear();
			std::getline( bndry_file, ln );
			// skip comment
			if ( ln.empty() || ln[ 0 ] == '#' )
				continue;

			ss.str( ln );
			ss.clear();
			ss >> v_id >> x >> y >> z;
			_pts[v_id] = point(x, y, z);
		}
		bndry_file.close();
		t_io.stop();
		cout << "time(I/O) -> read pts from .node: " << t_io.elapseMilli().count() << " ms" << endl;
		if ( v_id != n_bndry_vts - 1 )
		{
			cout << "Error: num of points read doesn't match the expected value!" << endl;
			return ImportErrCode::FAILURE;
		}

		return ImportErrCode::SUCCESS;
	}

	ImportErrCode voxelvoro::readVoroInfo( const char * _basename, VoroInfo& _voro, 
		bool _need_euler, const char* _sites_pos_file /*= nullptr*/ )
	{
		cout << "Reading voro info from tetgen files... " << endl;
		if ( !_sites_pos_file )
		{
			cout << "Error: cannot proceed because sites file is missing!" << endl;
			return ImportErrCode::FAILURE;
		}
		vector<point> sites_pos;
		cout << "Reading sites positions from file... " << endl;
		if ( readPts( _sites_pos_file, sites_pos ) != ImportErrCode::SUCCESS )
		{
			cout << "Error: failed to read sites position from file." << endl;
			return ImportErrCode::FAILURE;
		}
		cout << "Done: " << sites_pos.size() << " sites positions read. " << endl;
		// set closest sites position to voro vertices
		_voro.setSitesPositions( sites_pos );
		if ( _voro.loadFromTetgenFiles( _basename ) != true )
		{
			cout << "Error: couldn't load voro info from tetgen files." << endl;
			return ImportErrCode::FAILURE;
		}
		cout << "Done: voro info read." << endl;
		
		return ImportErrCode::SUCCESS;
	}
	ImportErrCode readMesh( const string& _filename, cellcomplex& _cc )
	{
		auto fname = fs::path( _filename );
		if ( fname.extension() == ".off" )
		{
			auto mesh_tmp = trimesh::TriMesh::read( _filename );
			if ( !mesh_tmp )
				return ImportErrCode::NOOPENINPUT;
			vector<uTriFace> trifaces;
			for ( auto f : mesh_tmp->faces )
				trifaces.emplace_back( f[ 0 ], f[ 1 ], f[ 2 ] );
			_cc = cellcomplex( mesh_tmp->vertices, {}, trifaces );
			delete mesh_tmp;
		}
		else if ( fname.extension() == ".ply" )
		{
			vector<float> dump;
			vector<point> vts;
			vector<ivec2> edges;
			vector<uTriFace> faces;
			readFromPLY( _filename.c_str(), vts, edges, faces, dump, dump, dump );
			_cc = cellcomplex( vts, edges, faces );
		}
		else
			return ImportErrCode::INVALID_FORMAT;
		return ImportErrCode::SUCCESS;
	}
	// TODO: add logic for reading V/E/F measures from file
	ImportErrCode readFromPLY( const char * _ply_filename, 
		vector<point>& _output_vts, vector<ivec2>& _output_edges, vector<uTriFace>& _output_tris, 
		vector<float>& _vts_msure, vector<float>& _edges_msure, vector<float>& _faces_msure )
	{
		std::map<string, PlyProperty> v_props_map;
		PlyProperty v_props[] = {// can add more properties (e.g. a scalar measure Vertex.s) 
			{ "x", Float32, Float32, offsetof( ply::Vertex, x ), 0, 0, 0, 0 },
			{ "y", Float32, Float32, offsetof( ply::Vertex, y ), 0, 0, 0, 0 },
			{ "z", Float32, Float32, offsetof( ply::Vertex, z ), 0, 0, 0, 0 } 
		};
		v_props_map[ "x" ] = v_props[ 0 ];
		v_props_map[ "y" ] = v_props[ 1 ];
		v_props_map[ "z" ] = v_props[ 2 ];
		map<string, PlyProperty> e_props_map;
		PlyProperty e_props[] = {
			{ "vertex1", Int32, Int32, offsetof( ply::Edge, v1 ), PLY_SCALAR, 0, 0, 0 },
			{ "vertex2", Int32, Int32, offsetof( ply::Edge, v2 ), PLY_SCALAR, 0, 0, 0 }
		};
		e_props_map[ "vertex1" ] = e_props[ 0 ];
		e_props_map[ "vertex2" ] = e_props[ 1 ];
		std::map<std::string, PlyProperty> f_props_map;
		PlyProperty f_props[] = {
			{ "vertex_indices", Int32, Int32, offsetof( ply::Face, verts ),
			PLY_LIST, Uint8, Uint8, offsetof( ply::Face,nvts ) }
		};
		f_props_map[ "vertex_indices" ] = f_props[ 0 ];
		vector<ply::Vertex> vts; // = *_skel_vts;
		vector<ply::Edge> edges; // = *_skel_edges;
		vector<ply::Face> faces; //*_skel_faces;

		ply::PLYreader ply_reader;
		auto err = ply_reader.read( _ply_filename, v_props_map, e_props_map, f_props_map, vts, edges, faces );
		if ( err != ply::SUCCESS )
		{
			cout << "Error: cannot read ply file: " << _ply_filename << endl;
			return ImportErrCode::NOOPENINPUT;
		}
		_output_vts.resize( vts.size() );
		_output_edges.resize( edges.size() );
		_output_tris.resize( faces.size() );
		// for now not reading any measure
		_vts_msure.clear(); _edges_msure.clear(); _faces_msure.clear();
		for ( auto i = 0; i < vts.size(); ++i )
		{
			const auto& v = vts[ i ];
			auto& u = _output_vts[ i ];
			u[ 0 ] = v.x;
			u[ 1 ] = v.y;
			u[ 2 ] = v.z;
		}
		for ( auto i = 0; i < edges.size(); ++i )
		{
			const auto& e = edges[ i ];
			auto& e1 = _output_edges[ i ];
			e1[ 0 ] = e.v1;
			e1[ 1 ] = e.v2;
		}
		for ( auto i = 0; i < faces.size(); ++i )
		{
			const auto& f = faces[ i ];
			auto& t = _output_tris[ i ];
			t[ 0 ] = f.verts[ 0 ];
			t[ 1 ] = f.verts[ 1 ];
			t[ 2 ] = f.verts[ 2 ];
		}
		return ImportErrCode::SUCCESS;
	}
	ImportErrCode readNumberList( const char * _nums_filename, vector<float>& _nums )
	{
		return ImportErrCode::SUCCESS;
	}
	ImportErrCode readMedialCurveInfo( 
		const char * _mc_geom_filename, 
		const char * _mc_msure_filename, 
		const char * _mc_order_filename, 
		const char* _skel_name,
		vector<point>& _mc_vts, vector<float>& _mc_msure, vector<int>& _mc_order,
		vector<point>* _skel_vts, vector<ivec2>* _skel_edges, vector<uTriFace>* _skel_faces )
	{
		std::ifstream mc_geom_os( _mc_geom_filename );
		if ( !mc_geom_os.is_open() )
		{
			std::cout << "Couldn't open MC geometry file: " << _mc_geom_filename << std::endl;
			return ImportErrCode::NOOPENINPUT;
		}
		int n_vts;
		mc_geom_os >> n_vts;
		_mc_vts.resize( n_vts );
		for ( auto& v : _mc_vts )
		{
			mc_geom_os >> v[ 0 ];
			mc_geom_os >> v[ 1 ];
			mc_geom_os >> v[ 2 ];
		}
		mc_geom_os.close();

		if ( _mc_msure_filename )
		{
			std::ifstream mc_msure_os( _mc_msure_filename );
			if ( !mc_msure_os.is_open() )
			{
				std::cout << "Couldn't open MC measure file: " << _mc_msure_filename << std::endl;
				return ImportErrCode::NOOPENINPUT;
			}
			int n_msure_entries;
			mc_msure_os >> n_msure_entries;
			assert( n_msure_entries == n_vts );
			_mc_msure.resize( n_msure_entries );
			for ( auto& s : _mc_msure )
				mc_msure_os >> s;
			mc_msure_os.close();
		}

		std::ifstream mc_order_os( _mc_order_filename );
		if ( !mc_order_os.is_open() )
		{
			std::cout << "Couldn't open MC vts' order file: " << _mc_order_filename << std::endl;
			return ImportErrCode::NOOPENINPUT;
		}
		int n_pairs;
		mc_order_os >> n_pairs;
		_mc_order.resize( n_pairs );
		for ( auto& o : _mc_order )
		{
			mc_order_os >> o >> o; // really just need the second number
		}
		mc_order_os.close();

		// optionally load skeleton vertices
		if ( _skel_name != nullptr && _skel_vts != nullptr )
		{
			_skel_vts->clear();
			_skel_edges->clear();
			_skel_faces->clear();
			vector<float> msure_dump;
			auto err = readFromPLY( _skel_name, *_skel_vts, *_skel_edges, *_skel_faces, msure_dump, msure_dump, msure_dump );
			if ( err != ImportErrCode::SUCCESS )
				return err;
		}
		
		return ImportErrCode::SUCCESS;
	}
}
