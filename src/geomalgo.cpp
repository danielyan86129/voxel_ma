#include "geomalgo.h"

#include <queue>
#include <unordered_set>
#include <unordered_map>

namespace util {
	using std::queue;
	using std::unordered_set;
	using std::unordered_map;

	void simpleTriangulate( const vector<int>& _poly, vector<uTriFace>& _tris )
	{
		_tris.clear();
		_tris.reserve( _poly.size() - 2 );
		int i = _poly.size() - 1, j = 0, k;
		bool inc_j = true;
		while ( i - 1 > j )
		{
			if ( inc_j )
			{
				k = j + 1;
				_tris.push_back( uTriFace( _poly[ i ], _poly[ j ], _poly[ k ] ) );
				j++;
			}
			else // decrease i
			{
				k = i - 1;
				_tris.push_back( uTriFace( _poly[ i ], _poly[ j ], _poly[ k ] ) );
				i--;
			}
			inc_j = !inc_j;
		}
	}
	void nonDegenTriangulate( const vector<point>& _poly, vector<uTriFace>& _tris )
	{

	}
	bool is_degenerate( const vector<point> _poly, float _eps )
	{
		// find a non-degenerate edge <u, v>
		int u, v;
		u = 0; v = 1;
		while ( trimesh::dist2( _poly[ u ], _poly[ v ] ) <= 0.0f )
		{
			if ( v == 0 )
				return true;
			u = v;
			v = ( v + 1 ) % _poly.size();
		}
		// if perpendicular distance from another vertex to uv is non-zero
		// then this is a non-degenerate poly
		auto uv_unit = _poly[ u ] - _poly[ v ];
		trimesh::normalize( uv_unit );
		for ( int w = 0; w < _poly.size(); ++w )
		{
			if ( w != u && w != v )
			{
				auto uw = _poly[ u ] - _poly[ w ];
				auto proj_len = uw.dot( uv_unit );
				auto h = std::sqrt( trimesh::len2( uw ) - trimesh::sqr( proj_len ) );
				if ( h > _eps )
					return false;
			}
		}
		return true;
	}
	int findNumConnComponents( const vector<point>& _vts, const vector<ivec2> _edges )
	{
		vector<vector<int>> vvadj( _vts.size(), {} );
		for ( auto e : _edges )
		{
			vvadj[ e[ 0 ] ].push_back( e[ 1 ] );
			vvadj[ e[ 1 ] ].push_back( e[ 0 ] );
		}
		vector<bool> visited( _vts.size(), false );
		return findNumConnComponents( vvadj, _vts.size(), visited );
	}
	int findNumConnComponents( 
		const vector<vector<int>>& _v_v_adj_tbl, 
		const vector<point>& _vts,
		vector<bool>& _visited )
	{
		int ncc = 0;
		ncc = findNumConnComponents( _v_v_adj_tbl, _vts.size(), _visited );

		return ncc;
	}
	int findNumConnComponents( const vector<vector<int>>& _v_v_adj_tbl, int _n_vts, vector<bool>& _visited )
	{
		int ncc = 0;
		int n_vts = _n_vts;

		queue<int> vts_q;
		for ( int vi = 0; vi < n_vts; ++vi )
		{
			if ( _visited[ vi ] )
				continue;
			// new component starts here
			ncc++;
			vts_q.push( vi );
			while ( !vts_q.empty() )
			{
				int cur_v = vts_q.front();
				vts_q.pop();
				if ( _visited[ cur_v ] )
					continue;
				_visited[ cur_v ] = true;
				// add its unvisited neighbors to q
				const auto& nbs = _v_v_adj_tbl[ cur_v ];
				for ( const auto& nb : nbs )
				{
					if ( !_visited[ nb ] )
						vts_q.push( nb );
				}
			}
		}

		return ncc;
	}
	void computeEulerChar( const vector<uTriFace>& _tris, eulerchar & _euler_struct )
	{
		vector<vector<int>> v_v_adj_tbl;
		unordered_set<ivec2, ivec2Hash> edges_unique;
		unordered_set<int> vts_unique_ids;
		int total_vts_size = 0;
		vector<bool> visited;
		ivec2 es[ 3 ];
		// collect unique edges and vts from faces info
		for ( const auto& f : _tris )
		{
			makeEdgesFromTri( f, es );
			for ( const auto& e : es )
			{
				edges_unique.insert( e );
				vts_unique_ids.insert( e[ 0 ] );
				vts_unique_ids.insert( e[ 1 ] );
			}
		}
		// set unreferenced vertices as visited already
		// (we don't want them to interfere with the major connected pieces)
		for ( auto vi : vts_unique_ids )
		{
			total_vts_size = std::max( total_vts_size, vi );
		}
		total_vts_size += 1;
		visited.resize( total_vts_size, false );
		for ( int vi = 0; vi < visited.size(); ++vi )
		{
			if (vts_unique_ids.count(vi) == 0 )
				visited[ vi ] = true;
		}
		// build v-v-adj tbl
		cout << "building v-v-adj table..." << endl;
		v_v_adj_tbl.resize( total_vts_size );
		for ( auto e_it = edges_unique.begin(); e_it != edges_unique.end(); ++e_it )
		{
			const auto& e = *e_it;
			v_v_adj_tbl[ e[ 0 ] ].push_back( e[ 1 ] );
			v_v_adj_tbl[ e[ 1 ] ].push_back( e[ 0 ] );
		}
		// num conn components
		cout << "computing num of connected components..." << endl;
		int ncc = findNumConnComponents( v_v_adj_tbl, vts_unique_ids.size(), visited );

		// fill-in euler struct
		int V = vts_unique_ids.size(); 
		int E = edges_unique.size();
		int F = _tris.size(); 
		int T = 0;
		int C = ncc;
		int euler = V - E + F - T;
		_euler_struct = eulerchar(V, E, F, T, C, euler);

		// clean up
		v_v_adj_tbl.clear();
		vts_unique_ids.clear();
		edges_unique.clear();
		visited.clear();
	}
	void computeEulerChar( 
		const vector<point>& _vts, const vector<ivec2>& _edges, const vector<uTriFace>& _tris, 
		eulerchar & _euler_struct )
	{
		std::unordered_set<ivec2, ivec2Hash> unique_edges;
		unique_edges.insert( _edges.begin(), _edges.end() );
		ivec2 tri_es[ 3 ];
		for ( auto& t : _tris )
		{
			makeEdgesFromTri( t, tri_es );
			unique_edges.insert( { tri_es[ 0 ], tri_es[ 1 ], tri_es[ 2 ] } );
		}
		int V = _vts.size();
		int E = unique_edges.size();
		int F = _tris.size();
		int T = 0;
		// # connect components
		// TODO: factor the logic to another function, e.g. connComps()
		vector<vector<int>> v_v_adj_tbl;
		unordered_set<int> vts_unique_ids;
		vector<bool> visited( V, false );
		v_v_adj_tbl.resize( V );
		for ( auto i = 0; i < V; ++i )
			vts_unique_ids.insert( i );
		for ( auto e_it = unique_edges.begin(); e_it != unique_edges.end(); ++e_it )
		{
			const auto& e = *e_it;
			v_v_adj_tbl[ e[ 0 ] ].push_back( e[ 1 ] );
			v_v_adj_tbl[ e[ 1 ] ].push_back( e[ 0 ] );
		}
		int C = findNumConnComponents( v_v_adj_tbl, vts_unique_ids.size(), visited );;
		int euler = V - E + F - T;
		_euler_struct = eulerchar( V, E, F, T, C, euler );
	}
	void compactify( const vector<int>& _v_l, vector<ivec2>& _e_l, vector<uTriFace>& _f_l )
	{
		// old-new-v-id-map
		unordered_map<int, int> old_new_map;
		for ( auto new_i = 0; new_i < _v_l.size(); ++new_i )
			old_new_map[ _v_l[ new_i ] ] = new_i;
		// rename e list
		for ( auto& e : _e_l )
		{
			e[ 0 ] = old_new_map[ e[ 0 ] ];
			e[ 1 ] = old_new_map[ e[ 1 ] ];
		}
		// rename f list
		for ( auto& f : _f_l )
		{
			f[ 0 ] = old_new_map[ f[ 0 ] ];
			f[ 1 ] = old_new_map[ f[ 1 ] ];
			f[ 2 ] = old_new_map[ f[ 2 ] ];
		}
	}

#define N_KEY_COLORS 9
	const TriColor jet_map[ N_KEY_COLORS ] = {
		TriColor( 0.0f, 0.0f, 0.6f ), // dark blue
		TriColor( 0.0f, 0.0f, 1.0f ), // blue
		TriColor( 0.0f, 0.5f, 1.0f ), // azure
		TriColor( 0.0f, 1.0f, 1.0f ), // cyan
		TriColor( 0.5f, 1.0f, 0.5f ), // light green
		TriColor( 1.0f, 1.0f, 0.0f ), // yellow
		TriColor( 1.0f, 0.5f, 0.0f ), // orange
		TriColor( 1.0f, 0.0f, 0.0f ), // red
		TriColor( 0.5f, 0.0f, 0.0f ) // dark red
	};

	TriColor GetColour( float v, float vmin, float vmax )
	{
		TriColor c( 1.0f, 1.0f, 1.0f ); // white
								
		// get the subrange that v falls into
		float range, dv;

		if ( v < vmin )
			v = vmin;
		if ( v > vmax )
			v = vmax;
		range = std::max(vmax - vmin, 0.00000000001f);
		dv = ( 1.0f / ( N_KEY_COLORS - 1 ) );

		float multiples = ( ( v - vmin ) / range ) / dv;
		int min_idx = multiples;
		const auto& min_c = jet_map[ min_idx ];
		const auto& max_c = jet_map[ min_idx + 1 ];
		c = trimesh::mix( min_c, max_c, ( multiples - min_idx ) );
		return( c );
	}

	bool intersect( const ivec3& _o1, const ivec3& _dim1, const ivec3& _o2, const ivec3& _dim2 )
	{
		auto min_1 = _o1; auto max_1 = _o1 + _dim1;
		auto min_2 = _o2; auto max_2 = _o2 + _dim2;
		bool x_overlap = !( min_1[ 0 ] > max_2[ 0 ] || min_2[ 0 ] > max_1[ 0 ] );
		bool y_overlap = !( min_1[ 1 ] > max_2[ 1 ] || min_2[ 1 ] > max_1[ 1 ] );
		bool z_overlap = !( min_1[ 2 ] > max_2[ 2 ] || min_2[ 2 ] > max_1[ 2 ] );

		return x_overlap && y_overlap && z_overlap;
	}
	
}// namespace geom