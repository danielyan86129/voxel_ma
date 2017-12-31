#ifndef VORO_INFO_H
#define VORO_INFO_H

#include <vector>
#include <memory>
#include <fstream>
#include <string>
#include <map>

//#include <volume.h> // Tao's volume
#include "Volume3DScalar.h"

#include "commondefs.h"
#include "spaceinfo.h"
#include "cellcomplex.h"
#include "measureforMA.h"
#include "ccthin.h"

namespace voxelvoro
{
	using std::vector;
	using std::shared_ptr;
	using std::iostream;
	using std::ifstream;
	using std::string;
	using std::map;

	enum class VoroRetCode
	{
		SUCCESS,
		FAILURE
	};

	/************************************************************************/
	/* Represents low-level info about a voro-diagram						*/
	/************************************************************************/
	class VoroInfo
	{
	public:
		VoroInfo();
		~VoroInfo();

		/*major interfaces*/
	public:
		//
		// load in voro info from related tetgen files with the given base name
		bool loadFromTetgenFiles( const char* _file_basename );
		//
		// init from given information
		void setInfo(
			const vector<point>& _vts, const vector<ivec2>& _edges, const vector<vector<int>> _faces,
			const vector<bool>& _vert_valid_tag,
			const vector<bool>& _is_finite_v,
			const vector<bool>& _is_finite_e,
			const vector<bool>& _is_finite_f
			);
		//
		// set closest sites to each voro vertex, 
		// and the radius function over the voro-diagram
		void setSitesPositions( const vector<point>& _sites_pos );
		// this version assumes sites for each face is correctly associated and only computes vert-radius 
		void computeInfoRelatedtoSites();
		// this version takes in the centroid for each cell and does 1-nn to associate sites with each face
		// then it also computes vert-radius
		void computeInfoRelatedtoSites( const vector<point>& _cell_cents );
		//
		// return positions of all sites
		const vector<point>& getSitesPosition() const;
		//
		// return whether radii over voro vertices are available
		bool radiiValid() const;
		//
		// return a list of radii defined over vertices of the voro-diagram
		// NOTE: radii only valid after they are initialized, e.g. by calling setClosestSites()
		const vector<float>& getRadii() const;
		//
		// return the contributing sites of a face
		ivec2 getSitesOfFace( int _fi ) const;
		//
		// tags each vert according to the voxel value of the given volume
		bool tagVtsUsingUniformVol( const shared_ptr<Volume3DScalar>& _vol );
		//
		// return the maximal side-length of the bbox of the inside part
		float getInsidePartSize() const;
		//
		// returns the list of valid vts
		void getValidVts( vector<int>& _vts_indices ) const;
		//
		// returns the list of valid edges (all vts must be valid)
		void getValidEdges( vector<int>& _edge_indices ) const;
		//
		// returns the list of valid faces (all vts must be valid)
		void getValidFaces( vector<int>& _face_indices ) const;
		//
		// returns the list of finite faces 
		void getFiniteFaces( vector<int>& _face_indices ) const;
		// 
		// returns the list of index for each face
		void getAllFaces( vector<int>& _face_indices ) const;
		// 
		// given indices, return actual geom representation (faces will be triangulated)
		void getVts( const vector<int>& _vts_indices, vector<point>& _vts ) const;
		void getEdges( const vector<int>& _edges_indices, vector<ivec2>& _edges ) const;
		void getFaces( const vector<int>& _faces_indices, vector<uTriFace>& tri_faces ) const;
		// 
		// get the voro-diagram i.e. vts, edges, and (triangulated)faces, with the invalid elements pruned. 
		void getPartInside( vector<int>& _vts_indices, vector<ivec2>& _edges, vector<uTriFace>& _tri_faces ) const;
		//
		// return inside part (v/e/f) together with the specified measure (on v/e/f)
		void extractInsideWithMeasure( MeasureForMA::meassuretype _mssure_tp,
			vector<point>& _output_vts, vector<ivec2>& _output_edges, vector<uTriFace>& _output_tris,
			vector<float>& _vts_msure, vector<float>& _edges_msure, vector<float>& _faces_msure ) const;
		//
		// 
		void generateMeasure( MeasureForMA::meassuretype _msure_type );
		//
		// return the geometry of this voro diagram
		inline const cellcomplex& geom() const
		{
			return m_geom;
		}
		//
		// merge those vertices close to each other
		// polygon faces and edges will be updated properly
		// @result current voro will be updated
		// @param _eps: within this range merging will occur
		// @param _only_inside: will only perform merging on inside elements & save IVD
		void mergeCloseVts( float _eps, bool _only_inside );
		//
		// return the euler characteristic for this voro
		void computeEulerChar( eulerchar& _euler_struct, bool _do_prune_first = true );
		//
		// compute the requested measure for vts/edges/faces
		void computeVertexMeasure( 
			MeasureForMA::meassuretype _mssure_tp, const vector<int>& _vts_indices, 
			vector<float>& _vts_msure ) const ;
		void computeEdgesMeasure( 
			MeasureForMA::meassuretype _mssure_tp, const vector<int>& _edges_indices, 
			vector<float>& _edges_msure ) const;
		void computeFacesMeasure( 
			MeasureForMA::meassuretype _mssure_tp, const vector<int>& _faces_indices, 
			vector<float>& _faces_msure ) const;
		// compute measure on edge as interpolation between those of end vts
		void computeEdgesMeasure(
			const vector<float>& _vts_msure, const vector<int>& _edges_indices,
			vector<float>& _edges_msure ) const;

		//
		// states computation / update
		inline bool computeEdgeValidity( int _ei ) const;
		inline bool computeEdgeValidity( const ivec2& _e ) const;
		inline bool computeFaceValidity( int _fi ) const;
		inline bool computeFaceValidity( const vector<int>& _f_vts ) const;
		inline bool computeFaceValidity( const uTriFace& _tri ) const;
		void recomputeEdgeFlags();
		void recomputeFaceFlags();

		//
		// states query interfaces
		inline bool isVertexValid( int _vi ) const;
		inline bool isEdgeValid( int _ei ) const;
		inline bool isFaceValid( int _fi ) const;
		inline bool isVertexFinite( int _vi ) const;
		inline bool isEdgeFinite( int _ei ) const;
		inline bool isFaceFinite( int _fi ) const;

		//
		// output the complex to a Mathematica-readable file
		void outputToMathematica( const char* _filename ) const;

	protected:
		/*helpers*/

		// 
		// invalidate geometry related states
		void invalidate_geometric_states();
		//
		// These only load in finite elements (i.e. finite edges, faces)
		void load_voro_vts( ifstream& _in_node );
		void load_voro_edges( ifstream& _in_edge, vector<ivec2>& _edges );
		void load_voro_faces( const vector<ivec2> &_edges, ifstream& _in_face );
		void infer_sites_from_cell_file( ifstream& _in_cells );
		//
		// trace a face boundary (a seq of edges) into a seq of vertices
		void trace_face( const vector<int>& _f_of_sides, const vector<ivec2>& _edges,
			vector<int>& _f_of_vts ) const;
		//
		// compute the multiplicity of each vertex, including the infinite vert (ref. by edge/face as -1)
		void compute_multiplicity_vts();

		/*
		protected data members
		*/
		// geometry stored in a cell-complex
		cellcomplex m_geom;
		// measure associated to each vert / edge / face
		vector<float> m_v_msure;
		vector<float> m_e_msure;
		vector<float> m_f_msure;
		// validity tag of each vert / edge/ face
		vector<bool> m_vts_valid;
		vector<bool> m_edge_valid;
		vector<bool> m_face_valid;
		// multiplicity of each vert (how many surface components touch at vi singularly)
		vector<int> m_mult_num_vts;
		// finiteness tag for each vert
		vector<bool> m_is_finite_v;
		// finiteness tag for each edge
		vector<bool> m_is_finite_e;
		// finiteness tag for each face
		vector<bool> m_is_finite_f;
		// the contributing sites (2) for each face
		vector<ivec2> m_face_sites;
		vector<point> m_site_positions;
		bool m_r_valid;// is radius function ready?
		bool m_face_sites_valid; // is sites for each face ready?
		// the radius function for the vertices
		vector<float> m_r_per_v;
		// the bbox of the finite part of the voro-diagram
		trimesh::box m_bbox;
	};
} // namespace voxelvoro

#include "voroinfo_imp.h"

#endif // #define VORO_INFO_H