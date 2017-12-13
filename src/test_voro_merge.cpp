#include <vector>
#include <memory>

#include "voroinfo.h"

using std::vector;
using std::shared_ptr;
using namespace voxelvoro;

vector<point> vts = {
	point( 0.0f, 0.0f, 0.0f ), // 0
	point( 1.0f, 1.0f, 0.0f ), // 1
	point( 1.0f, 1.0f, 0.0f ), // 2
	point( 0.0f, 2.0f, 0.0f ), // 3
	point( 2.0f, 1.0f, 0.0f ), // 4
	point( 2.0f, 1.0f, 0.0f ), // 5
	point( 3.0f, 2.0f, 0.0f ), // 6
	point( 3.0f, 0.0f, 0.0f ), // 7
	point( 0.0f, 2.0f, 0.0f ), // 8
	point( -1.0f, 1.0f, 0.0f ), // 9
};
// edges info.
vector<ivec2> edges = {
	// face 0
	ivec2( 0, 1 ), ivec2( 1, 2 ),ivec2( 2, 3 ),ivec2( 0, 3 ),
	// face 1
	ivec2( 2, 4 ), ivec2( 1, 4 ),
	// faec 2
	ivec2( 1, 5 ), ivec2( 4, 5 ),
	// face 3
	ivec2( 5, 7 ), ivec2( 6, 7 ),ivec2( 4, 6 ),
	// face 4
	ivec2( 0, 3 ), ivec2(3, 8), ivec2(0, 8),
	// face 5
	ivec2( 0, 8 ), ivec2( 8, 9 ), ivec2( 0, 9 )
};
vector<vector<int>> faces = {
	vector<int>( {0,1,2,3} ),
	vector<int>( { 1,2,4 } ),
	vector<int>( { 1,4,5 } ),
	vector<int>( { 4, 5, 7, 6 } ),
	vector<int>( {0, 3, 8} ),
	vector<int>( {0,8,9} )
};
vector<bool> vert_valid_tag = {
	true, true, true, true, true, true, true, true, true, true
};
vector<bool> is_finite_v = {
	true, true, true, true, true, true, true, true, true, true
};
// must indicate for each edge its finiteness
vector<bool> is_finite_e = {
	true, true, true, true, 
	true, true, 
	true, true,
	true, true, true,
	true, true, true,
	true, true, true
};
vector<bool> is_finite_f = {
	true, true, true, true, true, true
};

void main( void )
{
	VoroInfo voro;
	cout << "filling voro with given info ..." << endl;
	voro.setInfo( vts, edges, faces, vert_valid_tag, is_finite_v, is_finite_e, is_finite_f );
	cout << "Done: filling voro with given info" << endl;
	eulerchar euler;
	cout << "computing euler char ..." << endl;
	voro.computeEulerChar( euler );
	cout << "Done: computing euler char" << endl;
	euler.logToConsole( "pruned voro" );
	cout << "merging close vts ..." << endl;
	voro.mergeCloseVts( 0.00000001f, false );
	cout << "Done: merging close vts" << endl;
	cout << "computing euler char ..." << endl;
	voro.computeEulerChar( euler );
	cout << "Done: computing euler char" << endl;
	euler.logToConsole( "merged voro" );
}