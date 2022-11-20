#include <map>
#include <string>
#include <vector>

#include <ply/ply.h>

namespace ply
{
struct Vertex
{
    float x;
    float y;
    float z;
    unsigned char r, g, b;
    float s; // measure
};
struct Edge
{
    int v1;
    int v2;
    unsigned char r, g, b;
    float s; // measure
};
struct Face
{
    unsigned char nvts;
    int verts[3];
    unsigned char sites_l; // len of sites array
    float sites[6];        // 2*3 coords for two sites positions
    float s;               // measure
    unsigned char r, g, b; // pseudo-color of measure
};

enum class ErrCode
{
    SUCCESS,
    FAILURE
};

/// lower level reader / writer that knows underlying PLY parser format
/// from (ply lib)

class PLYReader
{
public:
    static ErrCode read(const char* _ply_filename,
                        const std::map<std::string, PlyProperty>& _v_props_map,
                        const std::map<std::string, PlyProperty>& _e_props_map,
                        const std::map<std::string, PlyProperty>& _f_props_map,
                        const std::vector<ply::Vertex>& _vts,
                        const std::vector<ply::Edge>& _edges,
                        const std::vector<ply::Face>& _faces);
};
class PLYWriter
{
public:
    static ErrCode write(const char* _ply_filename, bool _write_vts,
                         bool _write_edges, bool _write_faces,
                         const std::map<std::string, PlyProperty>& vert_props,
                         const std::map<std::string, PlyProperty>& _edge_props,
                         const std::map<std::string, PlyProperty>& _face_props,
                         const std::vector<Vertex>& output_vts,
                         const std::vector<Edge>& output_edges,
                         const std::vector<Face>& output_faces);
};
} // namespace ply