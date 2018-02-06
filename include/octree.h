#ifndef OCTREE_H
#define OCTREE_H

#include <string>
#include <vector>
//#include <volume.h> // Tao's volume class def.
#include "commondefs.h"
#include "Volume3DScalar.h"

using std::string;
using std::vector;

//////////////////////////////////////////////////////////////////////////
/// This octree represents a tree with the same layout as in a .sof file
//////////////////////////////////////////////////////////////////////////


/* the octree volume representation */
class OctreeVolume : public Volume3DScalar
{
public:
	struct OctreeNode;
	struct LeafNode;
	struct EmptyNode;
	struct IntNode;
public:
	// default constructor
	OctreeVolume();
	// constructor: read in octree from a .sof/.sog file
	OctreeVolume( const char* _file );
	// constructor: build octree from given volume
	OctreeVolume( const Volume3DScalar* _src_vol );
	// destructor
	~OctreeVolume();

	/*methods unique to OctreeVolume*/
	// write octree to file with known extension
	bool writeToFile( const char* _file ) const;
	// return the values for the requested node,
	// whose 8 corners represent 8 voxels
	unsigned char getNodeValues(int _i, int _j, int _k) const;
	// whether the given node is out side of the space represented by this octree
	bool outOfVolume( int _i, int _j, int _k ) const;
	// return a list of boundary voxels, 
	// i.e. those at the corners of leaf nodes with differing signs
	void getBoundaryVoxels( vector<ivec3>& _voxels ) const;
	// return the tightest bbox of the occupied volume
	trimesh::Box<3, int> getBoundingBox() const;

	/* static methods of OctreeVolume */
	// delete node recursively
	static void delnode( OctreeNode* _node )
	{
		if ( _node->type == 0 )
		{
			auto intnode = dynamic_cast<IntNode*>( _node );
			for ( auto i = 0; i < 8; ++i )
			{
				auto c = intnode->getChild( i );
				if ( c )
					delnode( c );
			}
		}
		delete _node;
	};

	/*implementing methods in parent class*/
	// Get data at a single voxel
	double getDataAt( int _x, int _y, int _z );
	double getDataAt( int x, int y, int z ) const;

protected: /*helpers*/
	//
	// helpers: read volume from a file of supported format (sof/sog)
	bool read_sof_file( const string& _sof_file );
	bool read_sog_file( const string& _sog_file );
	// write volume to a supported file format (sof/sog)
	bool write_sof_file( const string& _sof_file ) const;
	bool write_sog_file( const string& _sog_file ) const;
	// get the requested node's corners' values 
	// search within a node starting from _off of size _len
	// return whether the specified node is within the region or not
	bool get_node_values(
		int _i, int _j, int _k,
		OctreeNode* _node,
		ivec3 _off, int _len,
		unsigned char& _vals ) const;

	/*test*/
	// go through all nodes, test whether values returned by getDataAt()
	// match those stored represented by the node.
	bool test_getDataAt();
	// test if voxels values match the corners of a node
	bool test_node( OctreeNode* _node, ivec3 _off, int _len );

protected: /*static helpers*/
	// convert source volume into an octree node (covering a region with origin _p and extent _l)
	static OctreeNode* build( const ivec3& _p, int _l, const Volume3DScalar* _src_vol );

protected: /*data members*/
	// the file name
	std::string m_vol_filename;
	// max resolution of the octree (in terms of # nodes)
	// @Note the effective resolution for the voxels is this plus one
	int m_max_res;
	// the root of the octree
	OctreeNode* m_root;
	// the data with same layout as described in .sof file
	char* m_data;

	/*other octree related classes & interfaces*/
public:
	/*enumerate node types */
	enum NodeType
	{
		LEAF = 2, EMPTY = 1, INTERNAL = 0
	};
	/*Base type for a node*/
	struct OctreeNode
	{
		// 2: leaf, 1: empty, 0: internal
		unsigned char type;
		OctreeNode( unsigned char _type ) : type( _type ) {}
		virtual ~OctreeNode() {}
		OctreeVolume::NodeType getType() const { return ( OctreeVolume::NodeType )type; }
		virtual bool isHomo( unsigned char& _homo_val ) const { return true; }
	};
	/*a leaf node*/
	struct LeafNode : public OctreeNode
	{
		unsigned char value;
		LeafNode( unsigned char _val )
			: OctreeNode( NodeType::LEAF ), value( _val ) {}
		// if 8 corners have same values, return true and its value in homo_val
		// else return false
		bool isHomo( unsigned char& _homo_val ) const
		{
			_homo_val = value & 0x1; // any corner value would do
			return value == 0x00 || value == 0xff;
		}
		// get the value (0 | 1) at corner x-y-z
		int getCorner( int _x, int _y, int _z ) const
		{
			int corner_shift = getShift( _x, _y, _z );
			return ( ( value >> corner_shift ) & 1 ) == 0 ? 0 : 1;
		}
		// set the value of corner x-y-z to occval
		void setCorner( int _x, int _y, int _z, int _occval )
		{
			auto corner_shift = getShift( _x, _y, _z );
			value = value & ~( (unsigned char)( 1 ) << corner_shift ) | ( _occval << corner_shift );
		}
		// return an unchar to represent 8-voxels of this leafnode (for purpose of outputing to file)
		unsigned char getValues() const { return value; };
	protected:
		// return the shift-offset for getting the bit corresponding to the given corner
		int getShift( int _x, int _y, int _z ) const
		{
			int shift = ( 4 & ( _x << 2 ) ) | ( 2 & ( _y << 1 ) ) | _z;
			return shift;
		}
	};
	/*an empty node*/
	struct EmptyNode : public OctreeNode
	{
		unsigned char value;
		EmptyNode( unsigned char _val )
			: OctreeNode( NodeType::EMPTY ), value( _val ) {}
		bool isHomo( unsigned char& _homo_val ) const 
		{
			_homo_val = value;
			return true;
		}
		// return an unchar to represent 8-voxels of this leafnode (for purpose of outputing to file)
		unsigned char getValues() const { return value; };
	};
	/*an internal node*/
	struct IntNode : public OctreeNode
	{
		IntNode( ) : OctreeNode( NodeType::INTERNAL )
		{
			for ( auto c : children )
				c = nullptr;
		}
		bool isHomo( unsigned char& _homo_val ) const
		{
			return false;
		}
		// get the child corresponding to the given octant
		OctreeNode* getChild( int _i, int _j, int _k ) const
		{
			auto cidx = getMask( _i, _j, _k );
			return children[ cidx ];
		}
		OctreeNode* getChild( int _cidx ) const
		{
			return children[ _cidx ];
		}
		// set the child corresponding to the given octant
		void setChild( int _i, int _j, int _k, OctreeNode* _cp )
		{
			auto cidx = getMask( _i, _j, _k );
			children[ cidx ] = _cp;
		}
		void setChild( int _i, OctreeNode* _cp )
		{
			children[ _i ] = _cp;
		}
	protected:
		// children pointers
		OctreeNode* children[ 8 ];
		// return the index for getting the child corresponding to the given octant
		int getMask( int _i, int _j, int _k ) const
		{
			int childindex = ( 4 & ( _i << 2 ) ) | ( 2 & ( _j << 1 ) ) | _k;
			return childindex;
		}
	};

	/* a base class for *walker* from which the users should inherit 
	to make their own walker for specific purposes
	*/
	struct Walker
	{
		// process leaf node: either a LeafNode or EmptyNode
		virtual void leaf( OctreeVolume* _tree, OctreeNode* _node, const ivec3& _off, int _len ) {}
		// pre-process an internal node and returns whether pre-condition met or not
		virtual bool pre( OctreeVolume* _tree, OctreeNode* _node, const ivec3& _off, int _len )
		{
			return true;
		}
		// post-process an internal node and returns whether post-condition met or not
		virtual bool post( OctreeVolume* _tree, OctreeNode* _node, const ivec3& _off, int _len )
		{
			return true;
		}
	};

	// walk this tree with the given *walker* _w
	void walkTree( Walker* _w );
	// walk a subtree with the given *walker* _w
	void walkSubtree( Walker* _w, OctreeNode* _node, const ivec3& _off, int _len );
};

#endif // !OCTREE_H