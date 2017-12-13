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

/*Base type for a node*/
struct OctreeNode
{
	// 2: leaf, 1: empty, 0: internal
	unsigned char type;
	OctreeNode( unsigned char _type ) : type( _type ) {}
	virtual ~OctreeNode() {}
};
/*a leaf node*/
struct LeafNode : public OctreeNode
{
	unsigned char value;
	LeafNode( unsigned char _type, unsigned char _val ) 
		: OctreeNode( _type ), value( _val ) {}
	int getMask( int _i, int _j, int _k )
	{
		int cornerMask = ( 4 & ( _i << 2 ) ) | ( 2 & ( _j << 1 ) ) | _k;
		return cornerMask;
	}
};
/*an empty node*/
struct EmptyNode : public OctreeNode
{
	unsigned char value;
	EmptyNode( unsigned char _type, unsigned char _val )
		: OctreeNode( _type ), value( _val ) {}
};
/*an internal node*/
struct IntNode : public OctreeNode
{ 
	// children pointers
	OctreeNode* children[ 8 ];
	IntNode(unsigned char _type) : OctreeNode(_type)
	{
		for ( auto c : children )
			c = nullptr;
	}
	// get the child corresponding to the given octant
	OctreeNode* getChild( int _i, int _j, int _k )
	{
		auto cidx = ( 4 & ( _i << 2 ) ) | ( 2 & ( _j << 1 ) ) | _k;
		return children[ cidx ];
	}
	// set the child corresponding to the given octant
	void setChild( int _i, int _j, int _k, OctreeNode* _cp )
	{
		auto cidx = ( 4 & ( _i << 2 ) ) | ( 2 & ( _j << 1 ) ) | _k;
		children[ cidx ] = _cp;
	}
	void setChild( int _i, OctreeNode* _cp )
	{
		children[ _i ] = _cp;
	}
};

class OctreeVolume : public Volume3DScalar
{
public:
	// default constructor
	OctreeVolume();
	// constructor: read in octree from a .sof file
	OctreeVolume( const char* _sof_file );
	// destructor
	~OctreeVolume();

	/*methods unique to OctreeVolume*/
	// return the values for the requested node,
	// whose 8 corners represent 8 voxels
	unsigned char getNodeValues(int _i, int _j, int _k) const;
	// whether the given node is out side of the space represented by this octree
	bool outOfVolume( int _i, int _j, int _k ) const;
	// return a list of boundary voxels, 
	// i.e. those at the corners of leaf nodes with differing signs
	void getBoundaryVoxels( vector<ivec3>& _voxels ) const;

	/*implementing methods in parent class*/
	// Get data at a single voxel
	double getDataAt( int _x, int _y, int _z );
	double getDataAt( int x, int y, int z ) const;

private:
	/*helpers*/
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
	bool test_node( OctreeNode* _node );

private:
	// the file name
	std::string m_sof_file;
	// max resolution of the octree
	// @Note the effective resolution for the voxels is this plus one
	int m_max_res;
	// the root of the octree
	OctreeNode* m_root;
	// the data with same layout as described in .sof file
	char* m_data;
};

#endif // !OCTREE_H