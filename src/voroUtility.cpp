#include <filesystem>
#include <iostream>
using std::cout;
using std::endl;
#include <cstdio>
#include <cstring>
#include <string>
#include <map>
#include <gflags/gflags.h>

#include "importers.h"
#include "exporters.h"
#include "highlevelalgo.h"
#include "graphapp.h"

namespace fs = std::experimental::filesystem;
//
// Specify all flags referred by the main() entry later and their valid values.
// 
// define valid modes here
enum class ValidMode {
	MRC2MAT, VOL2MESH, VORO2MESH, VOL2MA, R, MRC2SOF, TREE, TOPO, TIME_IO
};
std::map<std::string, ValidMode> validmodes_map = {
	{"mrc2mat", ValidMode::MRC2MAT},
	{"vol2mesh", ValidMode::VOL2MESH},
	{"voro2mesh", ValidMode::VORO2MESH},
	{"vol2ma", ValidMode::VOL2MA},
	{"mrc2sof", ValidMode::MRC2SOF},
	{"r", ValidMode::R},
	{"t", ValidMode::TREE},
	{"topo", ValidMode::TOPO},
	{"timeio", ValidMode::TIME_IO}
};
//// flags/args for the program
//ValidMode mode; // current mode 
//bool need_euler; // need to compute euler-characteristic?
//// args for mode vol2mesh: 
//bool only_bndry_vts; // only output vts on the voxel boundary
//// args for mode voro2mesh:
//bool voro_out_2_mat; // output voro to a mathematica-friendly file
//bool voro_out_2_dotma; // output voro to .ma file

//
// define command-line flags
//

// define cmd flag 'md' for modes & its handler
DEFINE_string( md, "not-specified", "specify a mode for the utility to perform certain task. REQUIRED." );
static bool processMode( const char* _flagname, const std::string& _value )
{
	if ( validmodes_map.find( _value ) == validmodes_map.end() )
	{
		cout << "Invalid value for -md: " << _value << endl;
		return false;
	}
	return true;
}
DEFINE_validator( md, processMode );

// define cmd flag 'needEuler' 
DEFINE_bool( needEuler, false, "need to compute euler characteristic for the relevant geometry? OPTIONAL." );

// define cmd flag 'siteFile'
DEFINE_string( siteFile, "", "provide a .node file containing the contributing sites of a voro-diagram. OPTIONAL." );

// define cmd flag 'onlyBndryVts'
DEFINE_bool( onlyBndryVts, true, "only extract the vertices of the voxel boundary. OPTIONAL." );

// define cmd flag 'collapseEdges'
DEFINE_bool( collapseZeroLenEdges, false, "perform 0-length-edge collapse on the voro-diagram. OPTIONAL." );

// define cmd flag 'outToMat'
DEFINE_bool( outToMat, false, "output voro-diagram to a mathematica file. OPTIONAL" );

// define cmd flag 'outToDotma'
DEFINE_bool( outToDotma, false, "output voro-diagram to a .ma file (input format to QMAT). OPTIONAL." );

// cmd option: thinning threshold
DEFINE_string( tt, "0.08", "a list of comma-sep threshold used in thinning. OPTIONAL." );

// cmd option: mapping function from MC to voxel surface & related parameters
DEFINE_string( dofuncmap, "",
	"mapping specified function from MC to voxel surface: \n\
	[\"BT2\" | \"BT3\" | \"traveldist\" | \"seglabel\"] (Distabled by default). OPTIONAL." );
DEFINE_string( mcBase, "", "base name for a set of medial-curve related files. input to funcMapApp. OPTIONAL." );
DEFINE_double( smoothR, 0.1, "smoothing ratio for funcMapApp. OPTIONAL." );

// cmd options for -md=mrc2sof
DEFINE_string( mrc, "", "the input volume in format .mrc that's to be converted. REQUIRED." );
DEFINE_string( sof, "", "the output volume in format .sof the input will be converted to. OPTIONAL." );

// cmd options for -md=t
DEFINE_string( skm, "", "the input file containing graph and edge measures .skMsure. REQUIRED." );

/* cmd options for -md=topo */
DEFINE_string( meshToCheck, "",
	"The input file representing a mesh with open boundary (e.g. medial axis). REQUIRED."\
	"We want to detect whether there are any closed components - \"pockets\" on it." );
DEFINE_bool( isMan, true, "Is the meshToCheck a manifold? OPTIONAL." );
DEFINE_bool( loadIVD, true, "load only IVD from tetgen voro output. OPTIONAL." );

void printUsage()
{
	cout << google::ProgramUsage() << endl;
};

// split a string into a list of numbers
vector<float> split( const std::string& _s, char _d )
{
	vector<float> res;
	std::stringstream ss( _s );
	std::string item;
	while ( std::getline( ss, item, _d ) ) 
		res.push_back( std::stof( item ) );
	return res;
}

void main( int _argc, char * _argv[] )
{
	// set usage of the program
	google::SetUsageMessage( "\
		voroUtil -md=<a valid mode> args.\n \
		args depend on each mode. \n \
		" );

	// parse command line
	// flags (starting with -) will be removed. only the cmd-level arguments will remain.
	google::ParseCommandLineFlags( &_argc, &_argv, true );

	int n_remain_args = _argc - 1, cur_arg_idx = 1;

	// convert .mrc (volume) file to mathematica file
	if ( FLAGS_md == "mrc2mat" )
	{
		if ( n_remain_args < 2 )
		{
			cout << "wrong args: expecting .mrc/output file(s)" << endl;
			goto FAILURE;
		}
		else
		{
			shared_ptr<Volume3DScalar> vol;
			voxelvoro::readMRC( _argv[ cur_arg_idx ], vol );
			if ( voxelvoro::writeToMathematicaFromMRC( vol, _argv[ cur_arg_idx + 1 ] ) < 0 )
			{
				cout << "failed to convert files from " << _argv[ cur_arg_idx ] << " to " << _argv[ cur_arg_idx + 1 ] << endl;
			}
			cur_arg_idx += 2;
			goto SUCCESS;
		}
	}
	else if ( FLAGS_md == "vol2mesh" )
	{
		if ( n_remain_args < 2 )
		{
			cout << "wrong args: expecting a volume and an output mesh file name" << endl;
			goto FAILURE;
		}
		else
		{
			shared_ptr<Volume3DScalar> vol;
			auto retcode_vol = voxelvoro::readVolume( _argv[ cur_arg_idx ], vol );
			if ( retcode_vol == voxelvoro::ImportErrCode::SUCCESS )
			{
				cout << "Done: volume file reading." << endl;
			}
			else
			{
				cout << "Failed to read volume file. Quiting..." << endl;
				goto FAILURE;
			}

			auto retcode = voxelvoro::ExportErrCode::SUCCESS;
			if ( FLAGS_onlyBndryVts )
				voxelvoro::writeVolumeAsBoundaryPts( vol, _argv[ cur_arg_idx + 1 ] );
			else
				voxelvoro::writeVolumeAsBoundaryMesh( vol, _argv[ cur_arg_idx + 1 ], FLAGS_needEuler );
			if ( retcode != voxelvoro::ExportErrCode::SUCCESS )
			{
				cout << "failed to write boundary mesh from "
					<< _argv[ cur_arg_idx ] << " to " << _argv[ cur_arg_idx + 1 ] << endl;
				goto FAILURE;
			}
			cur_arg_idx += 2;
			goto SUCCESS;
		}
	}
	else if ( FLAGS_md == "r" )
	{
		// write radii field to file
		if ( n_remain_args < 3 )
		{
			cout << "wrong args: expecting medialaxis(.off), boundary(.node), and radii(output) files" << endl;
			goto FAILURE;
		}
		else
		{
			cout << "start estimating & writing radii field ..." << endl;

			int retcode = voxelvoro::estimateRadiiField(
				_argv[ cur_arg_idx ], _argv[ cur_arg_idx + 1 ], _argv[ cur_arg_idx + 2 ] );
			if ( retcode < 0 )
			{
				cout << "failed to estimate&write radii field to file." << endl;
				goto FAILURE;
			}
			cur_arg_idx += 3;

			cout << "Done estimating & writing radii field!" << endl;
			goto SUCCESS;
		}
	}
	else if ( FLAGS_md == "voro2mesh" )
	{
		if ( n_remain_args < 3 )
		{
			cout << "wrong args: expecting a volume file name, "
				<< "tetgen voro files basename, and an output mesh .off file name" << endl;
			goto FAILURE;
		}
		else
		{
			cout << "Reading a volume file: " << _argv[ cur_arg_idx ] << endl;
			shared_ptr<Volume3DScalar> vol;
			if ( voxelvoro::readVolume( _argv[ cur_arg_idx ], vol ) == voxelvoro::ImportErrCode::SUCCESS )
				cout << "Done: volume file reading." << endl;
			else
			{
				cout << "Error: couldn't read volume file." << endl;
				goto FAILURE;
			}

			auto basename = _argv[ cur_arg_idx + 1 ];
			cout << "Reading voroinfo ..." << endl;
			voxelvoro::VoroInfo voro;
			if (
				voxelvoro::readVoroInfo( basename, voro,
					FLAGS_needEuler, 
					FLAGS_siteFile == "" ? nullptr : FLAGS_siteFile.c_str(), 
					FLAGS_loadIVD ? vol : nullptr )
				== voxelvoro::ImportErrCode::SUCCESS )
				cout << "Done: voro info read." << endl;
			else
			{
				cout << "Error: couldn't read voro info." << endl;
				goto FAILURE;
			}

			cout << "Preprocessing voro... " << endl;
			if ( !voxelvoro::preprocessVoro( voro, vol, FLAGS_needEuler ) )
			{
				cout << "Error processing voro!" << endl;
				goto FAILURE;
			}
			cout << "Done: voro preprocessed." << endl;

			cout << "Writing inside voro to mesh file ..." << endl;
			auto voromesh_file = string( _argv[ cur_arg_idx + 2 ] );
			auto voromesh_filebase = voromesh_file.substr( 0, ( voromesh_file.find_last_of( '.' ) ) );
			if (
				voxelvoro::exportInsideVoroMesh( voro, vol, voromesh_file.c_str(),
					FLAGS_needEuler, FLAGS_collapseZeroLenEdges, 
					true/*inside only*/, false/*finite only*/, FLAGS_outToDotma,
					split(FLAGS_tt, ',') ) == voxelvoro::ExportErrCode::SUCCESS
				)
			{
				cout << "Done: inside voro info written." << endl;
			}
			else
			{
				cout << "Error: couldn't write inside voro info." << endl;
				goto FAILURE;
			}
			// optionally write out a radii file
			if ( voro.radiiValid() )
			{
				auto radii_filename = voromesh_filebase + ".r";
				cout << "writing radii to file " << radii_filename << endl;
				if ( voxelvoro::writeRadiiToFile( voro, radii_filename.c_str(), vol->getVoxToModelMat() )
					== voxelvoro::ExportErrCode::SUCCESS )
					cout << "Done: writing radii to file. " << endl;
				else
					cout << "Failed writing radii to file! " << endl;
			}

			cur_arg_idx += 3; n_remain_args -= 3;

			if ( FLAGS_dofuncmap != "" )
			{
				cout << endl << "***************************" << endl;
				cout << "Estimating & exporting scalar field on voxel surface..." << endl;
				auto mcgeom_filename = FLAGS_mcBase + ".mc";
				auto skel_filename = FLAGS_mcBase + "_skel.ply";
				if ( !std::experimental::filesystem::exists( skel_filename ) )
					skel_filename = "";
				cout << "skeleton file: " << skel_filename << endl;
				auto mcmsure_name = FLAGS_dofuncmap;
				auto mcmsure_filename = FLAGS_mcBase + "." + mcmsure_name + ".msure";
				// outputting segmentation needs help from another measure, e.g. "bt3"
				if ( mcmsure_name == "seglabel" || mcmsure_name == "length" )
					mcmsure_filename = FLAGS_mcBase + ".bt3.msure";
				else if ( mcmsure_name == "traveldist" )
					mcmsure_filename = FLAGS_mcBase + ".bt2.msure";
				cout << "measure file: " << mcmsure_filename << endl;
				auto mcorder_filename = FLAGS_mcBase + ".mcorder";
				auto output_filename = std::string( basename ) + "." + mcmsure_name + ".msure";
				auto retcode = voxelvoro::exportScalarFieldOnVoxelSurface(
					mcgeom_filename.c_str(),
					mcmsure_filename == "" ? nullptr : mcmsure_filename.c_str(),
					mcorder_filename.c_str(),
					skel_filename == "" ? nullptr : skel_filename.c_str(),
					output_filename.c_str(),
					mcmsure_name.c_str(), // name of the measure
					FLAGS_smoothR, // smoothing ratio (for smoothing measure)
					voro
				);
			}
			if ( FLAGS_outToMat )
			{
				// debug info required.
				auto matfilename = voromesh_filebase + "_voro_mat.txt";
				voro.outputToMathematica( matfilename.c_str() );
			}
			//if ( FLAGS_outToDotma )
			//{
			//	// output voro to a .ma file
			//	std::string vol_file = _argv[ 1 ];
			//	auto dotmafilename = vol_file.substr( 0, ( vol_file.find_last_of( '.' ) ) ) + ".ma";
			//	cout << "Done: writing voro-diagram to .ma file: " << dotmafilename << endl;
			//	if ( voxelvoro::writeInsideVoroToDotMA( voro, dotmafilename.c_str() )
			//		== voxelvoro::ExportErrCode::SUCCESS )
			//	{
			//		cout << "Done: writing voro-diagram to .ma file." << endl;
			//	}
			//	else
			//	{
			//		cout << "Error: failed to write voro-diagram to .ma file." << endl;
			//	}
			//}
			goto SUCCESS;
		}
	}
	else if ( FLAGS_md == "vol2ma" )
	{
		if ( n_remain_args < 2 )
		{
			cout << "wrong args: expecting <in: volume name> <out: MA meshname>" << endl;
			goto FAILURE;
		}
		else
		{
			cout << "Reading a volume file: " << _argv[ cur_arg_idx ] << endl;
			shared_ptr<Volume3DScalar> vol;
			if ( voxelvoro::readVolume( _argv[ cur_arg_idx ], vol ) == voxelvoro::ImportErrCode::SUCCESS )
				cout << "Done: volume file reading." << endl;
			else
			{
				cout << "Error: couldn't read volume file." << endl;
				goto FAILURE;
			}

			cout << "Computing voro ..." << endl;
			voxelvoro::VoroInfo voro;
			voxelvoro::computeVD( vol, voro );
			cout << "Done: voro computed." << endl;

			cout << "Preprocessing voro... " << endl;
			if ( !voxelvoro::preprocessVoro( voro, vol, FLAGS_needEuler ) )
			{
				cout << "Error processing voro!" << endl;
				goto FAILURE;
			}
			cout << "Done: voro preprocessed." << endl;

			cout << "Writing inside voro to mesh file ..." << endl;
			string vorocore_file = _argv[ cur_arg_idx + 1 ];
			auto vorocore_filebase = vorocore_file.substr( 0, ( vorocore_file.find_last_of( '.' ) ) );
			if (
				voxelvoro::exportInsideVoroMesh( voro, vol, vorocore_file.c_str(),
					FLAGS_needEuler, FLAGS_collapseZeroLenEdges, 
					true/*inside only*/, false/*finite only*/, FLAGS_outToDotma,
					split( FLAGS_tt, ',' ) ) == voxelvoro::ExportErrCode::SUCCESS
				)
			{
				cout << "Done: inside voro info written." << endl;
			}
			else
			{
				cout << "Error: couldn't write inside voro info." << endl;
				goto FAILURE;
			}
			// optionally write out a radii file
			if ( voro.radiiValid() )
			{
				auto radii_filename = vorocore_filebase + ".r";
				cout << "writing radii to file " << radii_filename << endl;
				if ( voxelvoro::writeRadiiToFile( voro, radii_filename.c_str(), vol->getVoxToModelMat() )
					== voxelvoro::ExportErrCode::SUCCESS )
					cout << "Done: writing radii to file. " << endl;
				else
					cout << "Failed writing radii to file! " << endl;
			}

			cur_arg_idx += 2; n_remain_args -= 2;

			if ( FLAGS_dofuncmap != "" )
			{
				cout << endl << "***************************" << endl;
				cout << "Estimating & exporting scalar field on voxel surface..." << endl;
				auto mcgeom_filename = FLAGS_mcBase + ".mc";
				auto skel_filename = FLAGS_mcBase + "_skel.ply";
				if ( !std::experimental::filesystem::exists( skel_filename ) )
					skel_filename = "";
				cout << "skeleton file: " << skel_filename << endl;
				auto mcmsure_name = FLAGS_dofuncmap;
				auto mcmsure_filename = FLAGS_mcBase + "." + mcmsure_name + ".msure";
				// outputting segmentation needs help from another measure, e.g. "bt3"
				if ( mcmsure_name == "seglabel" || mcmsure_name == "length" )
					mcmsure_filename = FLAGS_mcBase + ".bt3.msure";
				else if ( mcmsure_name == "traveldist" )
					mcmsure_filename = FLAGS_mcBase + ".bt2.msure";
				cout << "measure file: " << mcmsure_filename << endl;
				auto mcorder_filename = FLAGS_mcBase + ".mcorder";
				auto output_filename = std::string( vorocore_filebase ) + "." + mcmsure_name + ".msure";
				auto retcode = voxelvoro::exportScalarFieldOnVoxelSurface(
					mcgeom_filename.c_str(),
					mcmsure_filename == "" ? nullptr : mcmsure_filename.c_str(),
					mcorder_filename.c_str(),
					skel_filename == "" ? nullptr : skel_filename.c_str(),
					output_filename.c_str(),
					mcmsure_name.c_str(), // name of the measure
					FLAGS_smoothR, // smoothing ratio (for smoothing measure)
					voro
				);
			}
			if ( FLAGS_outToMat )
			{
				// debug info required.
				auto matfilename = vorocore_filebase + "_voro_mat.txt";
				voro.outputToMathematica( matfilename.c_str() );
			}
			//if ( FLAGS_outToDotma )
			//{
			//	// output voro to a .ma file
			//	std::string vol_file = _argv[ 1 ];
			//	auto dotmafilename = vol_file.substr( 0, ( vol_file.find_last_of( '.' ) ) ) + ".ma";
			//	cout << "Done: writing voro-diagram to .ma file: " << dotmafilename << endl;
			//	if ( voxelvoro::writeInsideVoroToDotMA( voro, dotmafilename.c_str() )
			//		== voxelvoro::ExportErrCode::SUCCESS )
			//	{
			//		cout << "Done: writing voro-diagram to .ma file." << endl;
			//	}
			//	else
			//	{
			//		cout << "Error: failed to write voro-diagram to .ma file." << endl;
			//	}
			//}
			goto SUCCESS;
		}
	}
	else if ( FLAGS_md == "mrc2sof" )
	{
		if ( !fs::exists( FLAGS_mrc ) )
		{
			cout << "Error: mrc file does not exists/or is not specified: " << FLAGS_mrc << endl;
			goto FAILURE;
		}
		auto mrc_path = fs::path( FLAGS_mrc );
		if ( mrc_path.extension() != ".mrc" )
		{
			cout << "Error: unsupported format for input volume: " << mrc_path.extension() << endl;
			goto FAILURE;
		}
		auto sof_path = fs::path( FLAGS_sof );
		if ( sof_path == "" )
		{
			sof_path = mrc_path;
			sof_path.replace_extension( ".sof" );
		}
		cout << "Converting...\n" << "input volume: " << mrc_path << "\n" << "output volume:\n" << sof_path << endl;
		if ( !voxelvoro::denseToSparse( mrc_path.string().c_str(), sof_path.string().c_str() ) )
			goto FAILURE;
		goto SUCCESS;
	}
	else if ( FLAGS_md == "t" )
	{
		vector<vector<float>> edge_msures;
		graphapp::WeightedGraph g;
		auto suc = graphapp::readGraph( FLAGS_skm, g, edge_msures );
		if ( suc )
			cout << "Done: graph built from file." << endl;
		else 
		{
			cout << "Failed to build graph from file." << endl;
			goto FAILURE;
		}
		vector<graphapp::NodeHandle> parent;
		graphapp::makeTreeFromGraph( g, edge_msures, graphapp::TreeMethod::ShortestPath, parent );
		cout << "Done: tree extracted from graph." << endl;

		string skm_tree_file( FLAGS_skm );
		auto pos = skm_tree_file.find( '.' );
		skm_tree_file.insert( pos, "_t" ); // indicating this is a tree
		skm_tree_file = fs::path( skm_tree_file ).replace_extension( ".ply" ).string();
		graphapp::exportTree( g, parent, skm_tree_file );
		cout << "Done: tree exported to file -> " << skm_tree_file << endl;

		goto SUCCESS;
	}
	else if ( FLAGS_md == "topo" )
	{
		// read in mesh
		cellcomplex cc;
		auto err = voxelvoro::readMesh( FLAGS_meshToCheck, cc, !FLAGS_isMan );
		if ( err != voxelvoro::ImportErrCode::SUCCESS )
		{
			cout << "Error: cannot read file: " << FLAGS_meshToCheck << endl;
			goto FAILURE;
		}
		// euler char
		eulerchar euler_char;
		cc.eulerChar( euler_char );
		cout << "euler char -> " << euler_char.euler << endl;
		cout << "connected component -> " << euler_char.C << endl;
		int n_closed = 0;
		if ( !FLAGS_isMan )
		{
			// remove open components. if anything left, the mesh has closed pockets.
			n_closed = voxelvoro::nClosedComponents( cc );
		}
		else
		{
			// # pockets equal to # conn. comp.
			n_closed = euler_char.C;
		}
		cout << "closed component -> " << n_closed << endl;
		goto SUCCESS;
	}
	else if ( FLAGS_md == "timeio" )
	{
		if ( n_remain_args < 1 )
		{
			cout << "wrong args: expecting tetgen voro files basename (i.e. w/o .v.*), " << endl;
			goto FAILURE;
		}
		else
		{
			auto basename = _argv[ cur_arg_idx ];
			cout << "Simply reading voro info to profile I/O time... " << endl;
			if ( voxelvoro::profileVoroInfoIO( basename ) != voxelvoro::ImportErrCode::SUCCESS )
			{
				cout << "Error: couldn't complete reading voro info!" << endl;
				goto FAILURE;
			}
			cout << "Done: I/O profiling." << endl;
		}
		goto SUCCESS;
	}
	else // unrecognized option
	{
		cout << "Error: option not recognized. This should never happen. " << endl;
		goto FAILURE;
	}

FAILURE:
	printUsage();
	exit( -1 );
SUCCESS:
	exit( 0 );
}