#include <set>
#include <unordered_map>
using std::set;
using std::unordered_map;
using std::pair;
#include <TriMesh.h>
#include "cellcomplex.h"

cellcomplex::cellcomplex()
{
	m_has_infinite_v = m_has_infinite_e = false;
}

cellcomplex::cellcomplex( const vector<point>& _vts, const vector<ivec2>& _edges, const vector<uTriFace>& _faces )
{
	for ( const auto& v : _vts )
		appendVert( v );
	for ( const auto& e : _edges )
		appendEdge( e );
	for ( const auto& t : _faces )
		appendFace( t );
	finalize();
}

cellcomplex::~cellcomplex()
{

}

void cellcomplex::appendFace( const vector<int>& _f_of_vts )
{
	int nsides = _f_of_vts.size();
	m_faceVtsN.push_back( nsides );
	auto& f_vrep = m_faces_v_rep[ nsides ];
	// save the index of the start of this face
	m_faceIndex.push_back( f_vrep.size() );
	// add to v-rep
	f_vrep.insert( f_vrep.end(), _f_of_vts.begin(), _f_of_vts.end() );
}

void cellcomplex::appendFace( const uTriFace & _f_tri )
{
	int nsides = 3;
	m_faceVtsN.push_back( nsides );
	auto& f_vrep = m_faces_v_rep[ nsides ];
	m_faceIndex.push_back( f_vrep.size() );
	f_vrep.insert( f_vrep.end(), _f_tri.begin(), _f_tri.end() );
}

void cellcomplex::appendFace(
	const vector<int>& _f_of_vts,
	const vector<int>& _f_of_sides )
{
	int nsides = _f_of_vts.size();
	m_faceVtsN.push_back( nsides );
	auto& f_vrep = m_faces_v_rep[ nsides ];
	auto& f_erep = m_faces_e_rep[ nsides ];
	// save the index of the start of this face
	m_faceIndex.push_back( f_vrep.size() );
	// add to v-rep & e-rep
	f_vrep.insert( f_vrep.end(), _f_of_vts.begin(), _f_of_vts.end() );
	f_erep.insert( f_erep.end(), _f_of_sides.begin(), _f_of_sides.end() );
}

size_t cellcomplex::createInfiniteVertex()
{
	if ( m_has_infinite_v )
		return m_inf_v_id;
	auto bbox = trimesh::box();
	for ( const auto& p : m_vts )
		bbox += p; 
	auto inf_vertex = bbox.center() + ( bbox.max - bbox.min );
	m_vts.push_back( inf_vertex );
	m_has_infinite_v = true;
	m_inf_v_id = m_vts.size() - 1;
	return m_inf_v_id;
}

size_t cellcomplex::createInfiniteEdge()
{
	if ( m_has_infinite_e )
		return m_inf_e_id;
	auto inf_v_id = createInfiniteVertex();
	m_edges.push_back( util::makeEdge( inf_v_id, inf_v_id ) );
	m_has_infinite_e = true;
	m_inf_e_id = m_edges.size() - 1;
	return m_inf_e_id;
}

size_t cellcomplex::createInfiniteVertex( const point & _p )
{
	m_vts.push_back( _p );
	m_inf_v_id = m_vts.size() - 1;
	m_has_infinite_v = true;
	return m_inf_v_id;
}

void cellcomplex::clear()
{
	clearVertList();
	clearEdgeList();
	clearFaceList();
	clear_adjacency();
}
void cellcomplex::clearFaceList()
{
	m_faceIndex.clear();
	m_faceVtsN.clear();
	m_faces_v_rep.clear();
	m_faces_e_rep.clear();
}
void cellcomplex::clearEdgeList()
{
	m_edges.clear();
}
void cellcomplex::clearVertList()
{
	m_vts.clear();
}

vector<int> cellcomplex::refCntPerVert() const
{
	vector<int> ref_cnt( numVts(), 0 );
	for ( auto i = 0; i < ref_cnt.size(); ++i )
	{
		ref_cnt[ i ] = cntNbEdgesofVert( i );
	}
	return ref_cnt;
}
vector<int> cellcomplex::refCntPerEdge() const
{
	vector<int> ref_cnt( numEdges(), 0 );
	for ( auto i = 0; i < ref_cnt.size(); ++i )
	{
		ref_cnt[ i ] = cntNbFacesofEdge( i );
	}
	return ref_cnt;
}

void cellcomplex::init_adjacency()
{
	m_nb_edges_for_v.assign( m_vts.size(), vector<int>() );
	m_nb_faces_for_e.assign( m_edges.size(), vector<int>() );
}

void cellcomplex::add_nb_edge( int _vi, int _ei )
{
	m_nb_edges_for_v[ _vi ].push_back( _ei );
}

void cellcomplex::add_nb_face( int _ei, int _nb_fi )
{
	m_nb_faces_for_e[ _ei ].push_back( _nb_fi );
}

void cellcomplex::clear_adjacency()
{
	m_nb_edges_for_v.clear();
	m_nb_faces_for_e.clear();
}


void cellcomplex::finalize()
{
	/* time stats */
	timer t_total;
	t_total.start();

	timer init_edges_t, extract_edges_t, save_edges_t;
	timer init_adj_t, build_adj_t;
	// check if the edges info has been finalized or not
	if ( m_faces_e_rep.empty() )
	{ // need to finalize 
		/*auto e_debug = geom::makeEdge( 1, 60 );
		vector<int> f_vts_debug = { 33,1,60,35 };*/
		// complete edge-idx rep. for faces
		unordered_map<ivec2, int, ivec2Hash> edge_index_map;
		vector<int> vts_of_f;
		vector<ivec2> edges_of_f;

		// init edge-index-map with existing edges before we add new edges
		// notice: order of edges will be kept as index
		init_edges_t.start();
		auto prev_num_of_edges = m_edges.size();
		for ( size_t ei = 0; ei < m_edges.size(); ++ei )
		{
			auto pre_size = edge_index_map.size();
			auto& e_idx = edge_index_map[ m_edges[ ei ] ];
			auto after_size = edge_index_map.size();
			if ( after_size > pre_size )
			{
				e_idx = pre_size;
			}
		}
		init_edges_t.stop();

		// extract edges used by faces
		extract_edges_t.start();
		for ( size_t fi = 0; fi < numFaces(); ++fi )
		{
			// we need to extract all edges of faces
			getFaceVRep( fi, vts_of_f );
			/*if ( vts_of_f == f_vts_debug )
				int stop = 1;*/
			util::traceFace( vts_of_f, edges_of_f );
			// save the edge-id rep. of the face into face-list
			for ( size_t ei = 0; ei < edges_of_f.size(); ++ei )
			{
				const auto& e = edges_of_f[ ei ];
				/*if ( e == e_debug )
					int stop = 1;*/
				int e_idx;
				auto find_it = edge_index_map.find( e );
				if ( find_it != edge_index_map.end() )
					e_idx = find_it->second;
				else
				{ // this is a new edge. remember its index.
					e_idx = edge_index_map.size();
					edge_index_map[ e ] = e_idx;
				}
				// actually save the edge-id
				m_faces_e_rep[ edges_of_f.size() ].push_back( e_idx );
			}
		}
		extract_edges_t.stop();

		// save each edge at its corresponding index in edge-list
		save_edges_t.start();
		m_edges.resize( edge_index_map.size() );
		for ( const auto& e_id_pair : edge_index_map )
		{
			/*if ( e_id_pair.first == e_debug )
				int stop = 1;*/
			m_edges[ e_id_pair.second ] = e_id_pair.first;
		}
		edge_index_map.clear();
		save_edges_t.stop();
	}

	// ready to build adjacency info
	init_adj_t.start();
	init_adjacency();
	init_adj_t.stop();

	// update adjacency info for vts
	build_adj_t.start();
	for ( size_t ei = 0; ei < numEdges(); ++ei )
	{
		const auto& e = getEdge( ei );
		add_nb_edge( e[ 0 ], ei );
		add_nb_edge( e[ 1 ], ei );
	}

	// update adjacency info for edges
	m_edge_ref_cnt.resize( m_edges.size(), 0 );
	vector<int> e_indices_of_f;
	for ( size_t fi = 0; fi < numFaces(); ++fi )
	{
		getFaceERep( fi, e_indices_of_f );
		for ( const auto& eid : e_indices_of_f )
		{
			m_edge_ref_cnt[ eid ] ++; // # ref by faces
			add_nb_face( eid, fi ); // nb face add to the edge's list
		}
	}
	build_adj_t.stop();
	t_total.stop();

	/*cout << "init edge struct: "
		<< init_edges_t.elapseMilli().count() << "ms" << endl;
	cout << "extract edges from faces: "
		<< extract_edges_t.elapseMilli().count() << "ms" << endl;
	cout << "store edges in order: "
		<< save_edges_t.elapseMilli().count() << "ms" << endl;
	cout << "init adj list: "
		<< init_adj_t.elapseMilli().count() << "ms" << endl;
	cout << "build adj list: "
		<< build_adj_t.elapseMilli().count() << "ms" << endl;*/
	cout << "time -> finalize CC:" << t_total.elapseMilli().count() << " ms" << endl;
}
