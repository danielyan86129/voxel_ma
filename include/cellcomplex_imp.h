#include "cellcomplex.h"
inline bool cellcomplex::isFinalized() const
{
	return m_is_finalized;
}
inline void cellcomplex::nbVts( unsigned _vi, vector<int>& _nbs ) const
{
	auto n_nbs = cntNbEdgesofVert( _vi );
	for ( auto i = 0; i < n_nbs; ++i )
	{
		auto e = getEdge( nbEdgeofVert( _vi, i ) );
		_nbs.push_back( e[ 0 ] == _vi ? e[ 1 ] : e[ 0 ] );
	}
}
inline int cellcomplex::cntNbEdgesofVert( unsigned _vi ) const
{
	return m_nb_edges_for_v[ _vi ].size();
}
inline int cellcomplex::cntNbFacesofEdge( unsigned _ei ) const
{
	return m_nb_faces_for_e[ _ei ].size();
}

inline int cellcomplex::nbEdgeofVert( unsigned _vi, int _ith_nb_e ) const
{
	return m_nb_edges_for_v[ _vi ][ _ith_nb_e ];
}

inline int cellcomplex::nbFaceofEdge( unsigned _ei, int _ith_nb_f ) const
{
	return m_nb_faces_for_e[ _ei ][ _ith_nb_f ];
}

inline vector<int> cellcomplex::nbEdgesofVert( unsigned _vi ) const
{
	return m_nb_edges_for_v[_vi];
}

inline vector<int> cellcomplex::nbFacesofEdge( unsigned _ei ) const
{
	return m_nb_faces_for_e[ _ei ];
}
