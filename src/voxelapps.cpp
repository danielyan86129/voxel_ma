#include "voxelapps.h"
#include "cellcomplex.h"
#include <queue>
//#include <TriMesh.h>
#include <TriMesh_algo.h>
#include <ANN/ANN.h>

namespace voxelvoro
{
	namespace apps
	{
		using std::queue;
		/************************
		** helpers
		*************************/
		//
		// return next neighbor of vi that's unlabled (label == 0), and lying on edge (not part of 2-face)
		int next_unlabeled_nb( const cellcomplex& _cc, const vector<int>& _skel_vts_label, int _vi )
		{
			int next = -1;
			// only care about isolated edge
			for ( auto ni = 0; ni < _cc.cntNbEdgesofVert( _vi ); ++ni )
			{
				if ( _cc.cntNbFacesofEdge( _cc.nbEdgeofVert( _vi, ni ) ) > 0 ) 
					continue;
				auto nb_e = _cc.getEdge( _cc.nbEdgeofVert( _vi, ni ) );
				auto other = ( nb_e[ 0 ] == _vi ? nb_e[ 1 ] : nb_e[ 0 ] );
				if ( _skel_vts_label[ other ] == 0 )
				{
					next = other;
					break;
				}
			}
			return next;
		};

		void segment( 
			const vector<point>& _skel_vts, 
			const vector<ivec2>& _skel_edges, 
			const vector<uTriFace>& _skel_faces, 
			vector<int>& _skel_vts_label )
		{
			_skel_vts_label.assign( _skel_vts.size(), 0 );
			int l = 1; // cur label
			cellcomplex cc( _skel_vts, _skel_edges, _skel_faces );
			cout << "Skel cc has vts/edges/faces = " << cc.numVts() << "/" << cc.numEdges() << "/" << cc.numFaces() << endl;

			// segmenting edges
			auto cnt_edge_conn_cmpnt = 0;
			for ( auto i = 0; i < cc.numVts(); ++i )
			{
				if ( cc.cntNbEdgesofVert( i ) != 2 )
				{
					//cout << "degree of v " << i << ": " << cc.cntNbEdgesofVert( i ) << endl;
					auto next = -1;
					// try to form a segment for each branch at this junction
					do {
						// skip if all nbs that are labeled already
						next = next_unlabeled_nb( cc, _skel_vts_label, i );
						if ( next == -1 )
							break;
						cout << "new seg found." << endl;
						cnt_edge_conn_cmpnt++;
						l++;
						// found the starting point of a new segment
						_skel_vts_label[ i ] = l; // assign cur vert label l
						_skel_vts_label[ next ] = l; // assign next vert label l too
						while ( cc.cntNbEdgesofVert( next ) == 2 ) // only continue if at interior point
						{
							next = next_unlabeled_nb( cc, _skel_vts_label, next );
							if ( next == -1 )
								break;
							_skel_vts_label[ next ] = l;
						}
					} while ( true );
				}
			}
			cout << "Done edge segmentation: " << cnt_edge_conn_cmpnt << " components." << endl;

			// segmenting faces (use flooding to identify components confined by non-manifold curves)
			vector<bool> visited_f( cc.numFaces(), false );
			queue<int> face_q;
			vector<int> f_erep; // edge-rep of a face
			vector<int> f_vrep; // vert-rep of a face
			auto cnt_face_conn_cmpnt = 0;
			for (auto i = 0; i < cc.numFaces(); ++i )
			{
				if ( visited_f[ i ] )
					continue;
				l++;
				cnt_face_conn_cmpnt++;
				face_q.push( i );
				while ( !face_q.empty() )
				{
					auto cur_f = face_q.front();
					face_q.pop();
					if ( visited_f[ cur_f ] )
						continue;
					visited_f[ cur_f ] = true;
					// mark vts of cur_f with cur label
					cc.getFaceVRep( cur_f, f_vrep );
					for ( auto vi : f_vrep )
						_skel_vts_label[ vi ] = l;
					// add unvisited neighbor faces to q
					cc.getFaceERep( cur_f, f_erep );
					for ( auto ei : f_erep )
					{
						if ( cc.cntNbFacesofEdge(ei) != 2 ) // stop propagating at non-manifold edge
							continue;
						for ( auto nb_fi : cc.nbFacesofEdge( ei ) )
						{
							if ( !visited_f[ nb_fi ] )
								face_q.push( nb_fi );
						}
					}
				}
			}
			cout << "Done face segmentation: " << cnt_face_conn_cmpnt << " components." << endl;

			// sanity check
			int zero_skel_vts = 0;
			for ( auto i = 0; i < _skel_vts_label.size(); ++i )
			{
				if ( _skel_vts_label[ i ] == 0 )
				{
					zero_skel_vts++;
					cout << "found 0-labeled skel v. degree = " << cc.cntNbEdgesofVert( i ) << endl;
				}
			}
			if ( zero_skel_vts )
				cout << "# skel vts having label 0: " << zero_skel_vts << endl;
		}

		//
		// associate each voxel site with the closest voro vert
		// and use the scalar of that voro vert
		void from_voxel_to_voro(
			const VoroInfo& _voro,
			vector<std::pair<int, float> >& _closest_voro_v_for_site )
		{
			auto& voxel_sites = _voro.getSitesPosition();
			std::unordered_set<int> nb_faces;

			_closest_voro_v_for_site.assign(
				voxel_sites.size(),
				{ -1, std::numeric_limits<float>::max() }
			);
			for ( auto i = 0; i < _voro.geom().numVts(); ++i )
			{
				if ( !_voro.isVertexValid( i ) )
					continue;
				const point& cur_voro_v = _voro.geom().getVert( i );
				// find the sites assoc.ed with this voro v
				nb_faces.clear();
				for ( const auto& ei : _voro.geom().nbEdgesofVert( i ) )
				{
					auto& fcs = _voro.geom().nbFacesofEdge( ei );
					nb_faces.insert( fcs.begin(), fcs.end() );
				}
				for ( auto fi : nb_faces )
				{
					auto sites = _voro.getSitesOfFace( fi );
					for ( auto j = 0; j < 2; ++j )
					{
						auto st = sites[ j ];
						auto& vi_d_pr = _closest_voro_v_for_site[ st ];
						auto d = ::trimesh::dist( voxel_sites[ st ], cur_voro_v );
						if ( d < vi_d_pr.second )
						{
							vi_d_pr.first = i;
							vi_d_pr.second = d;
						}
					}
				}
			}
		}

		void match_voro_with_medialcurve(
			const VoroInfo& _voro, 
			const vector<point>& _mc_vts,
			vector<int>& _voro_v_to_mc_v )
		{
			_voro_v_to_mc_v.assign( _voro.geom().numVts(), -1 );
			// insert MC vts to kd-tree
			ANNpointArray data_pts;
			ANNpoint query_p;
			ANNidxArray nn_idx = new ANNidx[ 1 ];
			auto dists = new ANNdist[ 1 ];
			ANNkd_tree* kdtree;
			query_p = annAllocPt( 3 );
			data_pts = annAllocPts( _mc_vts.size(), 3 );
			for ( auto i = 0; i < _mc_vts.size(); ++i )
			{
				const auto& v = _mc_vts[ i ];
				data_pts[ i ][ 0 ] = v[ 0 ];
				data_pts[ i ][ 1 ] = v[ 1 ];
				data_pts[ i ][ 2 ] = v[ 2 ];
			}
			kdtree = new ANNkd_tree( data_pts, _mc_vts.size(), 3 );
			std::cout << "Done: inserting MC vts to kdtree" << std::endl;
			// query each voro vert into kdtree to find corresponding MC vert
			for ( auto i = 0; i < _voro.geom().numVts(); ++i )
			{
				/*if ( i % 100000 )
				cout << "processed: " << i << endl;*/
				if ( !_voro.isVertexValid( i ) )
					continue;
				auto& v = _voro.geom().getVert( i );
				query_p[ 0 ] = v[ 0 ];
				query_p[ 1 ] = v[ 1 ];
				query_p[ 2 ] = v[ 2 ];
				kdtree->annkSearch( query_p, 1, nn_idx, dists );
				_voro_v_to_mc_v[ i ] = nn_idx[ 0 ];
			}
			std::cout << "Done: building correspondence between voro vts and MC vts." << std::endl;
			delete[] query_p;
			delete[] nn_idx;
			delete[] dists;
			annDeallocPts( data_pts );
			delete kdtree;
			annClose();
		}
		
		void from_mc_vert_to_mc_vert( 
			const vector<point>& _mc_vts, 
			const vector<float>& _mc_msure, MCMeasure _msure_type,
			const vector<int>& _mc_order,
			float _smooth_ratio,
			const vector<bool>* is_stable,
			vector<int>& _mc_v_to_mc_v,
			vector<float>& _smoothed_field )
		{
			_mc_v_to_mc_v.assign( _mc_vts.size(), 0 );
			_smoothed_field.assign( _mc_vts.size(), -1.0f );
			// determine smoothing threshold
			auto min_max_s = std::minmax_element( _mc_msure.begin(), _mc_msure.end() );
			float smooth_t = _smooth_ratio * ( *min_max_s.second - *min_max_s.first ) + *min_max_s.first;
			std::cout << "smoothing t = " << smooth_t << std::endl;

			for ( auto i = 0; i < _mc_vts.size(); ++i )
			{
				int cur = i, pre = i;
				float temp_s = _mc_msure[ i ];
				while ( _mc_order[ cur ] != cur )
				{
					temp_s += ::trimesh::dist( _mc_vts[ cur ], _mc_vts[ _mc_order[ cur ] ] );
					cur = _mc_order[ cur ];
					temp_s = std::max( temp_s, _mc_msure[ cur ] );
					if ( std::abs( temp_s - _mc_msure[ cur ] ) > smooth_t )
						break;
				}
				_mc_v_to_mc_v[ i ] = cur;
			}

			// finalize smoothed scalar field
			for ( auto i = 0; i < _mc_vts.size(); ++i )
			{
				int cur = i;
				// Optionally use the stable subset as guide, only stopping walking when it reaches a stable vert
				if ( is_stable )
					while ( !(*is_stable)[ cur ] ) cur = _mc_order[ cur ];
				//// just stop & use stable vert
				//smoothed_mc_field[ i ] = _mc_msure[ cur ];
				// or, we can allow some additional walking
				float d = _mc_msure[ cur ];
				if ( _msure_type == MCMeasure::BT3 | _msure_type == MCMeasure::BT2 )
				{
					d = _mc_msure[ _mc_v_to_mc_v[ _mc_v_to_mc_v[ cur ] ] ];
				}
				else if ( _msure_type == MCMeasure::RIDGE )
				{
					// scalar at i + d(i, i's destination)
					d = 0.0f;//_mc_msure[ i ];
					auto cur_i = cur, dest = _mc_v_to_mc_v[ cur ];
					while ( cur_i != dest )
					{
						d += ::trimesh::dist( _mc_vts[ cur_i ], _mc_vts[ _mc_order[ cur_i ] ] );
						cur_i = _mc_order[ cur_i ];
					}
				}
				else if ( _msure_type == MCMeasure::LENGTH )
				{
					d = _mc_msure[ i ];
					auto cur_i = i;
					while ( _mc_order[ cur_i ] != cur_i )
					{
						d += ::trimesh::dist( _mc_vts[ cur_i ], _mc_vts[ _mc_order[ cur_i ] ] );
						cur_i = _mc_order[ cur_i ];
					}
				}
				_smoothed_field[ i ] = d;
			}
			std::cout << "Done: smoothing." << std::endl;
		}

		//
		// tag a subset of medial-curve points as stable using skeleton vts
		// return tag in _is_stable 
		// & _mc_to_skel, that maps each stable mc vert to its corresponding vert (-1 means no stable match).
		void tag_stable_subset_with_skel(
			const VoroInfo& _voro,
			const vector<point>& _mc_vts,
			const vector<point>& _skel_vts,
			vector<bool>& _is_stable, // each mc vert: has a corresponding skel vert or not
			vector<int>& _mc_to_skel // mapping: mc vert -> skel vert 
		)
		{
			ANNpointArray data_pts;
			ANNpoint query_p;
			ANNidxArray nn_idx;
			auto dists = new ANNdist[ 1 ];
			ANNkd_tree* kdtree;
			query_p = annAllocPt( 3 );
			data_pts = annAllocPts( _mc_vts.size(), 3 );
			for ( auto i = 0; i < _mc_vts.size(); ++i )
			{
				const auto& v = _mc_vts[ i ];
				data_pts[ i ][ 0 ] = v[ 0 ];
				data_pts[ i ][ 1 ] = v[ 1 ];
				data_pts[ i ][ 2 ] = v[ 2 ];
			}
			kdtree = new ANNkd_tree( data_pts, _mc_vts.size(), 3 );

			// actual tagging, use skeleton vts as queries
			_mc_to_skel.assign( _mc_vts.size(), -1 );
			auto eps = _voro.getInsidePartSize() * 0.00000001f;
			for ( auto i = 0; i < _skel_vts.size(); ++i )
			{
				//if ( i % 100 ) cout << "processed skeleton points: " << i << endl; // debug
				auto sv = _skel_vts[ i ];
				query_p[ 0 ] = sv[ 0 ];
				query_p[ 1 ] = sv[ 1 ];
				query_p[ 2 ] = sv[ 2 ];
				nn_idx = new ANNidx[ 1 ];
				kdtree->annkSearch( query_p, 1, nn_idx, dists );
				auto k_expected = kdtree->annkFRSearch( query_p, dists[ 0 ] + eps, 0 );
				if ( k_expected == 0 ) // dbg
					cout << "Warning: skel vertex " << i << " has labeled 0 MC vts as stable!" << endl;
				delete[] nn_idx;
				nn_idx = new ANNidx[ k_expected ];
				k_expected = kdtree->annkFRSearch( query_p, dists[ 0 ] + eps, k_expected, nn_idx );
				for ( auto j = 0; j < k_expected; ++j )
				{
					_is_stable[ nn_idx[ j ] ] = true;
					_mc_to_skel[ nn_idx[ j ] ] = i;
					//// record the label if want to project "segment label" onto surface
					//if ( _msure_type == MCMeasure::SegmentLabel )
					//	smoothed_mc_field[ nn_idx[ j ] ] = skel_vts_label[ i ];
				}
				delete[] nn_idx;
			}
			delete[] query_p;
			delete[] dists;
			annDeallocPts( data_pts );
			delete kdtree;
			annClose();
			std::cout << "Done: labeling stable subset of MC." << std::endl;
		}

		void assignScalarToSites(
			const VoroInfo& _voro,
			const vector<point>& _mc_vts, const vector<float>& _mc_msure, const vector<int>& _mc_order,
			MCMeasure _msure_type,
			float _smooth_ratio,
			const CellComplexThinning* const _ccthin, /*optional*/
			const vector<point>* const _skel_vts,
			const vector<ivec2>* const _skel_edges,
			const vector<uTriFace>* const _skel_faces,
			vector<float>& _site_scalar )
		{
			// store the smoothed scalar field on MC here
			vector<float> smoothed_mc_field( _mc_vts.size(), -1.0f );

			// stable subset of MC vts indicator (by default all stationary points are stable)
			// will be updated & used in smoothing if skeleton vertices are given, additionally
			vector<bool> is_stable( _mc_vts.size(), false );
			for ( auto i = 0; i < _mc_vts.size(); ++i )
				is_stable[ i ] = ( _mc_order[ i ] == i );
			vector<int> mc_to_skel; // map each stable mc vert to its corresponding vert. -1 means no stable match.

			if ( _skel_vts )
			{
				tag_stable_subset_with_skel( _voro, _mc_vts, *_skel_vts, is_stable, mc_to_skel );
			}

			/* -------smooth the scalar field------- */
			vector<int> mc_to_mc;
			from_mc_vert_to_mc_vert( 
				_mc_vts, _mc_msure, _msure_type, _mc_order, _smooth_ratio, 
				_skel_vts ? &is_stable : nullptr, 
				mc_to_mc,
				smoothed_mc_field );
			// debug: output smoothed scalar field on MC.
			std::ofstream dbg_smoothed_os( "dbg_mc_smoothed.msure" );
			dbg_smoothed_os << smoothed_mc_field.size() << std::endl;
			for ( auto s : smoothed_mc_field )
			{
				dbg_smoothed_os << s << std::endl;
			}
			dbg_smoothed_os.close();

			/* -------build correspondence between voro's vts and MC vts------- */
			vector<int> voro_v_to_mc_v( _voro.geom().numVts(), -1 );
			match_voro_with_medialcurve( _voro, _mc_vts, voro_v_to_mc_v );
			
			/* -------assign smoothed scalar to contributing sites------- */
			auto& voxel_sites = _voro.getSitesPosition();
			_site_scalar.assign( voxel_sites.size(), -1.0f );
			vector<std::pair<int, float> > closest_voro_v_for_site;
			from_voxel_to_voro( _voro, closest_voro_v_for_site );
			cout << "Done: from voxel sites to voro vts." << endl;
			
			if ( _msure_type == MCMeasure::SegmentLabel )
			{
				/* Optional: segment skeleton structure if provided */
				vector<int> skel_vts_label( _skel_vts->size(), 0 );
				cout << "skel vts label length = " << skel_vts_label.size() << endl;
				voxelvoro::apps::segment( *_skel_vts, 
					_skel_edges ? *_skel_edges : vector<ivec2>(), 
					_skel_faces ? *_skel_faces : vector<uTriFace>(), 
					skel_vts_label );
				// sanity check: how many segments MC will be partitioned into
				std::unordered_set<int> labels;
				for ( auto i : mc_to_skel )
					if ( i >= 0 )
						labels.insert( skel_vts_label[ i ] );
				cout << "Done checking. # segments on MC: " << labels.size() << endl;

				for ( auto i = 0; i < voxel_sites.size(); ++i )
				{
					int dest_mc_v = voro_v_to_mc_v[ closest_voro_v_for_site[ i ].first ];
					while ( mc_to_skel[ dest_mc_v ] < 0 ) dest_mc_v = _mc_order[ dest_mc_v ];
					if ( _smooth_ratio > 0 )
						dest_mc_v = mc_to_mc[ dest_mc_v ];
					_site_scalar[ i ] = skel_vts_label[ mc_to_skel[ dest_mc_v ] ];
				}
			}
			else // other measures
			{
				for ( auto i = 0; i < voxel_sites.size(); ++i )
				{
					_site_scalar[ i ] = smoothed_mc_field[ voro_v_to_mc_v[ closest_voro_v_for_site[ i ].first ] ];
				}
			}
		}
	}
}