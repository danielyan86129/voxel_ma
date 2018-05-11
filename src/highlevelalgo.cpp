#include "highlevelalgo.h"

#include <chrono>
#include <sstream>

#include <TriMesh.h>
#include <TriMesh_algo.h>
#include <plyall.h>

#include "surfacing.h"
#include "octree.h"
#include "voxelapps.h"

using trimesh::TriMesh;

namespace voxelvoro
{
	void my_tetrahedralize( char * switches, tetgenio * _in, tetgenio* out, 
		tetgenmesh& m, tetgenio * addin, tetgenio * bgmin )
	{
		tetgenbehavior b;
		if ( !b.parse_commandline( switches ) )
		{
			terminatetetgen( nullptr, 10 );
		}
		my_tetrahedralize( &b, _in, out, m, addin, bgmin );
	}
	void my_tetrahedralize( tetgenbehavior * _b, tetgenio * _in, tetgenio* _out, 
		tetgenmesh& _m, tetgenio * _addin, tetgenio * _bgmin )
	{
		clock_t tv[ 12 ], ts[ 5 ]; // Timing informations (defined in time.h)
		REAL cps = (REAL)CLOCKS_PER_SEC;

		tv[ 0 ] = clock();

		_m.b = _b;
		_m.in = _in;
		_m.addin = _addin;

		if ( _b->metric && _bgmin && ( _bgmin->numberofpoints > 0 ) ) {
			_m.bgm = new tetgenmesh(); // Create an empty background mesh.
			_m.bgm->b = _b;
			_m.bgm->in = _bgmin;
		}

		_m.initializepools();
		_m.transfernodes();

		exactinit( _b->verbose, _b->noexact, _b->nostaticfilter,
			_m.xmax - _m.xmin, _m.ymax - _m.ymin, _m.zmax - _m.zmin );

		tv[ 1 ] = clock();

		if ( _b->refine ) { // -r
			_m.reconstructmesh();
		}
		else { // -p
			_m.incrementaldelaunay( ts[ 0 ] );
		}

		tv[ 2 ] = clock();

		if ( !_b->quiet ) {
			if ( _b->refine ) {
				printf( "Mesh reconstruction seconds:  %g\n", ( (REAL)( tv[ 2 ] - tv[ 1 ] ) ) / cps );
			}
			else {
				printf( "Delaunay seconds:  %g\n", ( (REAL)( tv[ 2 ] - tv[ 1 ] ) ) / cps );
				if ( _b->verbose ) {
					printf( "  Point sorting seconds:  %g\n", ( (REAL)( ts[ 0 ] - tv[ 1 ] ) ) / cps );
				}
			}
		}

		if ( _b->plc && !_b->refine ) { // -p
			_m.meshsurface();

			ts[ 0 ] = clock();

			if ( !_b->quiet ) {
				printf( "Surface mesh seconds:  %g\n", ( (REAL)( ts[ 0 ] - tv[ 2 ] ) ) / cps );
			}

			if ( _b->diagnose ) { // -d
				_m.detectinterfaces();

				ts[ 1 ] = clock();

				if ( !_b->quiet ) {
					printf( "Self-intersection seconds:  %g\n", ( (REAL)( ts[ 1 ] - ts[ 0 ] ) ) / cps );
				}

				// Only output when self-intersecting faces exist.
				if ( _m.subfaces->items > 0l ) {
					_m.outnodes( _out );
					_m.outsubfaces( _out );
				}

				return;
			}
		}

		tv[ 3 ] = clock();

		if ( ( _b->metric ) && ( _m.bgm != NULL ) ) { // -m
			_m.bgm->initializepools();
			_m.bgm->transfernodes();
			_m.bgm->reconstructmesh();

			ts[ 0 ] = clock();

			if ( !_b->quiet ) {
				printf( "Background mesh reconstruct seconds:  %g\n",
					( (REAL)( ts[ 0 ] - tv[ 3 ] ) ) / cps );
			}

			if ( _b->metric ) { // -m
				_m.interpolatemeshsize();

				ts[ 1 ] = clock();

				if ( !_b->quiet ) {
					printf( "Size interpolating seconds:  %g\n", ( (REAL)( ts[ 1 ] - ts[ 0 ] ) ) / cps );
				}
			}
		}

		tv[ 4 ] = clock();

		if ( _b->plc && !_b->refine ) { // -p
			if ( _b->nobisect ) { // -Y
				_m.recoverboundary( ts[ 0 ] );
			}
			else {
				_m.constraineddelaunay( ts[ 0 ] );
			}

			ts[ 1 ] = clock();

			if ( !_b->quiet ) {
				if ( _b->nobisect ) {
					printf( "Boundary recovery " );
				}
				else {
					printf( "Constrained Delaunay " );
				}
				printf( "seconds:  %g\n", ( (REAL)( ts[ 1 ] - tv[ 4 ] ) ) / cps );
				if ( _b->verbose ) {
					printf( "  Segment recovery seconds:  %g\n", ( (REAL)( ts[ 0 ] - tv[ 4 ] ) ) / cps );
					printf( "  Facet recovery seconds:  %g\n", ( (REAL)( ts[ 1 ] - ts[ 0 ] ) ) / cps );
				}
			}

			_m.carveholes();

			ts[ 2 ] = clock();

			if ( !_b->quiet ) {
				printf( "Exterior tets removal seconds:  %g\n", ( (REAL)( ts[ 2 ] - ts[ 1 ] ) ) / cps );
			}

			if ( _b->nobisect ) { // -Y
				if ( _m.subvertstack->objects > 0l ) {
					_m.suppresssteinerpoints();

					ts[ 3 ] = clock();

					if ( !_b->quiet ) {
						printf( "Steiner suppression seconds:  %g\n",
							( (REAL)( ts[ 3 ] - ts[ 2 ] ) ) / cps );
					}
				}
			}
		}

		tv[ 5 ] = clock();

		if ( _b->coarsen ) { // -R
			_m.meshcoarsening();
		}

		tv[ 6 ] = clock();

		if ( !_b->quiet ) {
			if ( _b->coarsen ) {
				printf( "Mesh coarsening seconds:  %g\n", ( (REAL)( tv[ 6 ] - tv[ 5 ] ) ) / cps );
			}
		}

		if ( ( _b->plc && _b->nobisect ) || _b->coarsen ) {
			_m.recoverdelaunay();
		}

		tv[ 7 ] = clock();

		if ( !_b->quiet ) {
			if ( ( _b->plc && _b->nobisect ) || _b->coarsen ) {
				printf( "Delaunay recovery seconds:  %g\n", ( (REAL)( tv[ 7 ] - tv[ 6 ] ) ) / cps );
			}
		}

		if ( ( _b->plc || _b->refine ) && _b->insertaddpoints ) { // -i
			if ( ( _addin != NULL ) && ( _addin->numberofpoints > 0 ) ) {
				_m.insertconstrainedpoints( _addin );
			}
		}

		tv[ 8 ] = clock();

		if ( !_b->quiet ) {
			if ( ( _b->plc || _b->refine ) && _b->insertaddpoints ) { // -i
				if ( ( _addin != NULL ) && ( _addin->numberofpoints > 0 ) ) {
					printf( "Constrained points seconds:  %g\n", ( (REAL)( tv[ 8 ] - tv[ 7 ] ) ) / cps );
				}
			}
		}

		if ( _b->quality ) {
			_m.delaunayrefinement();
		}

		tv[ 9 ] = clock();

		if ( !_b->quiet ) {
			if ( _b->quality ) {
				printf( "Refinement seconds:  %g\n", ( (REAL)( tv[ 9 ] - tv[ 8 ] ) ) / cps );
			}
		}

		if ( ( _b->plc || _b->refine ) && ( _b->optlevel > 0 ) ) {
			_m.optimizemesh();
		}

		tv[ 10 ] = clock();

		if ( !_b->quiet ) {
			if ( ( _b->plc || _b->refine ) && ( _b->optlevel > 0 ) ) {
				printf( "Optimization seconds:  %g\n", ( (REAL)( tv[ 10 ] - tv[ 9 ] ) ) / cps );
			}
		}

		if ( !_b->nojettison && ( ( _m.dupverts > 0 ) || ( _m.unuverts > 0 )
			|| ( _b->refine && ( _in->numberofcorners == 10 ) ) ) ) {
			_m.jettisonnodes();
		}

		if ( ( _b->order == 2 ) && !_b->convex ) {
			_m.highorder();
		}

		if ( !_b->quiet ) {
			printf( "\n" );
		}

		if ( _out != (tetgenio *)NULL ) {
			_out->firstnumber = _in->firstnumber;
			_out->mesh_dim = _in->mesh_dim;
		}

		if ( _b->nonodewritten || _b->noiterationnum ) {
			if ( !_b->quiet ) {
				printf( "NOT writing a .node file.\n" );
			}
		}
		else {
			_m.outnodes( _out );
		}

		if ( _b->noelewritten ) {
			if ( !_b->quiet ) {
				printf( "NOT writing an .ele file.\n" );
			}
			_m.indexelements();
		}
		else {
			if ( _m.tetrahedrons->items > 0l ) {
				_m.outelements( _out );
			}
		}

		if ( _b->nofacewritten ) {
			if ( !_b->quiet ) {
				printf( "NOT writing an .face file.\n" );
			}
		}
		else {
			if ( _b->facesout ) {
				if ( _m.tetrahedrons->items > 0l ) {
					_m.outfaces( _out );  // Output all faces.
				}
			}
			else {
				if ( _b->plc || _b->refine ) {
					if ( _m.subfaces->items > 0l ) {
						_m.outsubfaces( _out ); // Output boundary faces.
					}
				}
				else {
					if ( _m.tetrahedrons->items > 0l ) {
						_m.outhullfaces( _out ); // Output convex hull faces.
					}
				}
			}
		}


		if ( _b->nofacewritten ) {
			if ( !_b->quiet ) {
				printf( "NOT writing an .edge file.\n" );
			}
		}
		else {
			if ( _b->edgesout ) { // -e
				_m.outedges( _out ); // output all mesh edges. 
			}
			else {
				if ( _b->plc || _b->refine ) {
					_m.outsubsegments( _out ); // output subsegments.
				}
			}
		}

		if ( ( _b->plc || _b->refine ) && _b->metric ) { // -m
			_m.outmetrics( _out );
		}

		if ( !_out && _b->plc &&
			( ( _b->object == tetgenbehavior::OFF ) ||
			( _b->object == tetgenbehavior::PLY ) ||
				( _b->object == tetgenbehavior::STL ) ) ) {
			_m.outsmesh( _b->outfilename );
		}

		if ( !_out && _b->meditview ) {
			_m.outmesh2medit( _b->outfilename );
		}


		if ( !_out && _b->vtkview ) {
			_m.outmesh2vtk( _b->outfilename );
		}

		if ( _b->neighout ) {
			_m.outneighbors( _out );
		}

		if ( _b->voroout ) {
			_m.outvoronoi( _out );
		}


		tv[ 11 ] = clock();

		if ( !_b->quiet ) {
			printf( "\nOutput seconds:  %g\n", ( (REAL)( tv[ 11 ] - tv[ 10 ] ) ) / cps );
			printf( "Total running seconds:  %g\n", ( (REAL)( tv[ 11 ] - tv[ 0 ] ) ) / cps );
		}

		if ( _b->docheck ) {
			_m.checkmesh( 0 );
			if ( _b->plc || _b->refine ) {
				_m.checkshells();
				_m.checksegments();
			}
			if ( _b->docheck > 1 ) {
				_m.checkdelaunay();
			}
		}

		if ( !_b->quiet ) {
			_m.statistics();
		}
	}
	void pts2tetgen( const vector<point> _P, tetgenio& _tetio )
	{
		if ( _tetio.numberofpoints != _P.size() )
		{
			if ( _tetio.pointlist )
				delete [] _tetio.point2tetlist;
			_tetio.pointlist = new REAL[ 3 * _P.size() ];
		}
		for ( auto i = 0; i < _P.size(); ++i )
		{
			auto p = _P[ i ];
			_tetio.pointlist[ i * 3 + 0 ] = p[ 0 ];
			_tetio.pointlist[ i * 3 + 1 ] = p[ 1 ];
			_tetio.pointlist[ i * 3 + 2 ] = p[ 2 ];
		}
		_tetio.numberofpoints = _P.size();
	}
	void computeVD( const shared_ptr<Volume3DScalar>& _vol, VoroInfo & _voro )
	{
		timer t;
		vector<point> sites;
		Surfacer surf;

		t.start();
		cout << "extracting sites from voxel boundary..." << endl;
		surf.extractBoundaryVts( _vol, sites );
		printf( "Done: extracted sites %d \n", sites.size() );
		tetgenio tet_io;
		pts2tetgen( sites, tet_io );
		t.stop();
		cout << "time -> preparing for VD: " << t.elapseMilli().count() << endl;
		computeVD( tet_io, _voro, _vol );
	}
	void computeVD( tetgenio & _tet_in, VoroInfo & _voro, shared_ptr<Volume3DScalar> _vol )
	{
		char* tet_args = "NEFv"; // don't need tets, but want voronoi
		tetgenio tet_out;
		/*tetgenmesh tet_msh;
		my_tetrahedralize( tet_args, &_tet_in, &tet_out, tet_msh );*/
		tetrahedralize( tet_args, &_tet_in, &tet_out );
		
		// set sites to voro
		vector<point> sites;
		for ( auto i = 0; i < _tet_in.numberofpoints; ++i )
		{
			sites.emplace_back( 
				_tet_in.pointlist[ i * 3 ],
				_tet_in.pointlist[ i * 3 + 1 ],
				_tet_in.pointlist[ i * 3 + 2 ] );
		}
		_voro.setSitesPositions( sites );
		sites.clear();
		_voro.loadFromTetgenFiles( tet_out, _vol );
	}
	bool preprocessVoro( VoroInfo & _voro, const shared_ptr<Volume3DScalar>& _vol, bool _need_euler )
	{
		timer t_tagging, t_merging, t_compute_msure;

		if ( !_voro.onlyHasInside() )
		{
			t_tagging.start();
			if ( !_voro.tagVtsUsingUniformVol( _vol ) )
			{
				cout << "Error: couldn't tag voro vts." << endl;
				return false;
			}
			t_tagging.stop();
		}

		// compute measures for voro elements
		t_compute_msure.start();
		//std::cout << "Computing measures ..." << std::endl;
		_voro.generateMeasure( MeasureForMA::LAMBDA );
		//std::cout << "Done: computing measures" << std::endl;
		t_compute_msure.stop();

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
		
		//cout << "time -> compute measure: " << t_compute_msure.elapseMilli().count() << " ms" << endl;
		cout << "time -> VD in/out tagging: " << t_tagging.elapseMilli().count() << " ms" << endl;
		cout << "time(FULLSTAGE) -> IVD-merge: " << t_merging.elapseMilli().count() << " ms" << endl;
		if ( _need_euler )
		{	// check euler numbers after merging identical vts
			_voro.computeEulerChar( euler_struct );
			euler_struct.logToConsole(
				( string( "after merging vts within " ) + std::to_string( eps ) ).c_str() );
		}
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
		const VoroInfo& _voro )
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
		int _full_or_thin, /*0: thin, 1: full, 2: both*/
		bool _write_attrib,
		bool _out_qmat,
		const vector<float>& _tt/*thinning threshold(s)*/ )
	{
		timer t_thin, t_write_IVD, t_write_thinIVD;

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

		auto mesh_file = std::string( _mesh_file );
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
		else if ( mesh_file.find( ".ply" ) != std::string::npos )
		{
			string mesh_path_wo_ext = mesh_file.substr( 0, mesh_file.rfind( '.' ) );
			vector<point> vts;
			vector<ivec2> edges;
			vector<uTriFace> tri_faces;
			vector<int> from_fi;
			vector<float> vts_msure, edges_msure, all_tri_msure;
			_voro.extractInsideWithMeasure( MeasureForMA::LAMBDA,
				vts, edges, tri_faces, from_fi,
				vts_msure, edges_msure, all_tri_msure );
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

				// output full VC to file
				t_write_IVD.start();
				if ( _full_or_thin == 1 || _full_or_thin == 2 ) // writing full VC?
				{
					cout << "writing full voxel core (could take a while) to: " << _mesh_file << endl;
					if (_write_attrib )
					{
						// prepare sites for all faces
						vector<ivec2> triface_site_ids;
						for ( auto i = 0; i < tri_faces.size(); ++i )
							triface_site_ids.push_back( _voro.getSitesOfFace( from_fi[ i ] ) );
						// write them along with full VC
						writeToPLY( _mesh_file, vts_trans, edges, tri_faces,
							vts_msure, edges_msure, all_tri_msure,
							true, &_voro.getSitesPosition(), &triface_site_ids );
					}
					else
					{
						writeToPLY( _mesh_file, vts_trans, edges, tri_faces,
							vts_msure, edges_msure, all_tri_msure );
					}
					auto radii_filename = mesh_path_wo_ext + ".r";
					cout << "writing radii to file " << radii_filename << endl;
					if ( voxelvoro::writeRadiiToFile( _voro, radii_filename.c_str(), _vol->getVoxToModelMat() )
						== voxelvoro::ExportErrCode::SUCCESS )
						cout << "Done: writing radii to file. " << endl;
					else
						cout << "Failed: writing radii to file! " << endl;
				}
				t_write_IVD.stop();
				auto min_max_msure = std::minmax_element( all_tri_msure.begin(), all_tri_msure.end() );
				cout << "range of measure: " << "[" << *min_max_msure.first << "," << *min_max_msure.second << "]" << endl;

				if ( _full_or_thin == 1 ) // not writing thinned VC?
					goto END_PLY_BRANCH;

				// qmat file related vars
				unordered_set<ivec2, ivec2Hash> unique_edges;

				std::stringstream ss;
				float eps = _voro.getInsidePartSize() * 0.000001f;
				CellComplexThinning ccthin;
				cellcomplex cc( vts_trans, edges, tri_faces );
				ccthin.setup( &cc );
				ccthin.assignElementValues( vts_msure, edges_msure, all_tri_msure );
				ccthin.preprocess();
				vector<int> remain_vts_ids;
				for ( double t : _tt )
				{
					//if ( util::is_equal( t, 0.0, (double)eps ) && _full_or_thin )
					//	continue; // skip this as full VC already written above!

					// start writing this pruned VC
					if ( t >= 0.0 )
					{
						// output pruned mesh
						t_thin.start();

						float thresh = _voro.getInsidePartSize() * t;
						ccthin.prune( thresh, thresh, false );
						ccthin.remainingCC( vts_trans, edges, tri_faces, nullptr );

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

						//* output the thinned result to .ply file
						std::string thin_file = _mesh_file;
						ss.str( "" ); 
						ss << t;
						auto end_str = std::string( "_thinned" ) + ss.str();
						thin_file.insert( thin_file.find( ".ply" ), end_str );
						t_write_thinIVD.start();
						if ( _write_attrib )
						{
							// write attrib for each face: measure & 2 site positions...
							vector<int> remain_tri_ids;
							ccthin.getRemainedFaces( remain_tri_ids );
							vector<ivec2> triface_site_ids;
							vector<float> tri_msure_thin;
							for ( auto ti : remain_tri_ids )
							{
								triface_site_ids.push_back( _voro.getSitesOfFace( from_fi[ ti ] ) ); // sites for remain tri faces
								tri_msure_thin.push_back( all_tri_msure[ ti ] ); // msure for remain tri faces
							}
							cout << "writing remaining cc to file (with face attribs): " << thin_file << endl;
							vts_msure.clear(); edges_msure.clear();
							writeToPLY( thin_file.c_str(), vts_trans, edges, tri_faces,
								vector<float>(), vector<float>(), tri_msure_thin,
								true, &( _voro.getSitesPosition() ), &triface_site_ids );
							cout << "Done." << endl;
						}
						else
						{
							cout << "writing remaining cc to file: " << thin_file << endl;
							writeToPLY( thin_file.c_str(), vts_trans, edges, tri_faces,
								vector<float>(), vector<float>(), vector<float>(),
								false );
							cout << "Done." << endl;
						}
						t_write_thinIVD.stop();
						// output .r file
						auto radii_filename = mesh_path_wo_ext + end_str + ".r";
						remain_vts_ids.clear();
						ccthin.getRemainedVts( remain_vts_ids );
						voxelvoro::writeRadiiToFile( _voro, radii_filename.c_str(), _vol->getVoxToModelMat(), &remain_vts_ids );

						//* output the thinned result to qmat .ma file
						if ( _out_qmat )
						{
							// get radii
							vector<int> vts_ids;
							ccthin.getRemainedVts( vts_ids );
							vector<float> radii;
							const auto& R = _voro.getRadii();
							auto orig_trans = _vol->getVoxToModelMat() * point( 0, 0, 0 );
							for ( auto i : vts_ids )
							{
								// transform r
								float r = trimesh::len( ( _vol->getVoxToModelMat() * point( R[ i ], 0, 0 ) ) - orig_trans );
								radii.push_back( r );
							}
							// assemble all edges
							unique_edges.clear();
							ivec2 es_f[ 3 ];
							for ( const auto& f : tri_faces )
							{
								util::makeEdgesFromTri( f, es_f );
								unique_edges.insert( es_f[ 0 ] );
								unique_edges.insert( es_f[ 1 ] );
								unique_edges.insert( es_f[ 2 ] );
							}
							edges.insert( edges.end(), unique_edges.begin(), unique_edges.end() );
							unique_edges.clear();
							// write to qmat
							auto dotmafilename = thin_file.substr( 0, ( thin_file.find_last_of( '.' ) ) ) + ".ma";
							auto dotma_suc = writeToDotma( dotmafilename, vts_trans, edges, tri_faces, radii );
							if ( dotma_suc == voxelvoro::ExportErrCode::SUCCESS )
							{
								cout << "Done: writing thinned voro-diagram to .ma file: " << dotmafilename << endl;
							}
							else
							{
								cout << "Error: failed to write thinned voro-diagram to .ma file." << endl;
							}
						}
					} // pruning for one threshold
				} // pruning for all thresholds in _tt list
			} // local scope for goto to skip
		END_PLY_BRANCH:
			cout << "time(THIN) -> thinning: " << t_thin.elapseMilli().count() << " ms" << endl;
			cout << "time(I/O) -> write IVD: " << t_write_IVD.elapseMilli().count() << " ms" << endl;
			cout << "time(I/O) -> write thinned IVD: " << t_write_thinIVD.elapseMilli().count() << " ms" << endl;
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