#include <vector>
#include <memory>

#include "edgecollapse.h"

using std::vector;
using std::shared_ptr;
using namespace voxelvoro;

//point vts[] = {
//	point( 0.0f, 0.0f, 0.0f ),
//	point( 1.0f, 0.0f, 0.0f ),
//	point( 1.0f, 0.0f, 0.0f ),
//	point( 1.0f, 0.0f, 0.0f )
//};
//uTriFace tris[] = {
//	uTriFace( 0, 1, 2 ),
//	uTriFace( 0, 1, 3 )
//};
point vts[] = {
	point( 0.0f, 0.0f, 0.0f ),
	point( 1.0f, 0.0f, 0.0f ),
	point( 1.0f, 0.0f, 0.0f ),
	point( 0.0f, 0.0f, 0.0f )
};
uTriFace tris[] = {
	uTriFace( 0, 1, 2 ),
	uTriFace( 1, 2, 3 )
};
int main( void )
{
	// build vts and faces lists
	shared_ptr<vector<point>> vts_ptr = std::make_shared<vector<point>>();
	shared_ptr<vector<uTriFace>> tris_ptr = std::make_shared<vector<uTriFace>>();
	for ( const auto& v : vts )
	{
		vts_ptr->push_back( v );
	}
	for ( const auto& f : tris )
	{
		tris_ptr->push_back( f );
	}
	shared_ptr<EdgeCollapser> ec = std::make_shared<TopoPreservEdgeCollapser>( vts_ptr, tris_ptr );
	ec->preProcess();
	ec->collapse();

	return 0;
}