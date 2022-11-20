#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

#include <voxelcore/plyall.h>

namespace ply
{
template <typename PtrType>
void* cast_to_nonconst_void_ptr(PtrType* x)
{
    return const_cast<void*>(static_cast<const void*>(x));
}

ErrCode PLYReader::read(const char* _ply_filename,
                        const std::map<std::string, PlyProperty>& _v_props_map,
                        const std::map<std::string, PlyProperty>& _e_props_map,
                        const std::map<std::string, PlyProperty>& _f_props_map,
                        const std::vector<Vertex>& _vts,
                        const std::vector<Edge>& _edges,
                        const std::vector<Face>& _faces)
{
    float version = 0.1;
    FILE* fp = fopen(_ply_filename, "r");
    if (fp == nullptr)
    {
        std::cout << "Failed to open file for reading " << _ply_filename
                  << std::endl;
        return ErrCode::FAILURE;
    }

    PlyFile* ply_file = read_ply(fp);
    for (int i = 0; i < ply_file->num_elem_types; ++i)
    {
        int elem_count = 0;
        char* elem_name = setup_element_read_ply(ply_file, i, &elem_count);

        if (equal_strings("vertex", elem_name))
        {
            std::vector<Vertex> vts;
            setup_property_ply(
                ply_file, const_cast<PlyProperty*>(&(_v_props_map.at("x"))));
            setup_property_ply(
                ply_file, const_cast<PlyProperty*>(&(_v_props_map.at("y"))));
            setup_property_ply(
                ply_file, const_cast<PlyProperty*>(&(_v_props_map.at("z"))));

            // read in all vertices
            for (int j = 0; j < elem_count; ++j)
            {
                get_element_ply(ply_file, cast_to_nonconst_void_ptr(&vts[j]));
            }
        }
        else if (equal_strings("edge", elem_name))
        {
            std::vector<Edge> edges;
            setup_property_ply(ply_file, const_cast<PlyProperty*>(
                                             &(_e_props_map.at("vertex1"))));
            setup_property_ply(ply_file, const_cast<PlyProperty*>(
                                             &(_e_props_map.at("vertex2"))));

            // read in all edges
            for (int j = 0; j < elem_count; ++j)
            {
                get_element_ply(ply_file, cast_to_nonconst_void_ptr(&edges[j]));
            }
        }
        else if (equal_strings("face", elem_name))
        {
            std::vector<Edge> faces;
            setup_property_ply(
                ply_file,
                const_cast<PlyProperty*>(&(_f_props_map.at("vertex_indices"))));

            // read in all faces
            for (int j = 0; j < elem_count; ++j)
            {
                get_element_ply(ply_file, cast_to_nonconst_void_ptr(&faces[j]));
            }
        }
        else
        {
            // skip this unknown data by "reading"
            std::cout << "====> Skipping unknown element type: " << elem_name
                      << std::endl;
            get_other_element_ply(ply_file);
        }
    }
    // that's it
    close_ply(ply_file);
    free_ply(ply_file);

    return ErrCode::SUCCESS;
}

ErrCode PLYWriter::write(const char* _ply_filename, bool _write_vts,
                         bool _write_edges, bool _write_faces,
                         const std::map<std::string, PlyProperty>& _vert_props,
                         const std::map<std::string, PlyProperty>& _edge_props,
                         const std::map<std::string, PlyProperty>& _face_props,
                         const std::vector<Vertex>& _output_vts,
                         const std::vector<Edge>& _output_edges,
                         const std::vector<Face>& _output_faces)
{
    FILE* fp = fopen(_ply_filename, "w");
    if (fp == nullptr)
    {
        std::cout << "Failed to open file for writing " << _ply_filename
                  << std::endl;
        return ErrCode::FAILURE;
    }

    std::vector<char*> elem_names = {"vertex", "edge", "face"};
    PlyFile* ply_file =
        write_ply(fp, elem_names.size(), elem_names.data(), PLY_ASCII);

    // setup vertex props
    describe_element_ply(ply_file, "vertex", _output_vts.size());
    describe_property_ply(ply_file,
                          const_cast<PlyProperty*>(&(_vert_props.at("x"))));
    describe_property_ply(ply_file,
                          const_cast<PlyProperty*>(&(_vert_props.at("y"))));
    describe_property_ply(ply_file,
                          const_cast<PlyProperty*>(&(_vert_props.at("z"))));

    // setup edge props
    describe_element_ply(ply_file, "edge", _output_edges.size());
    describe_property_ply(
        ply_file, const_cast<PlyProperty*>(&(_edge_props.at("vertex1"))));
    describe_property_ply(
        ply_file, const_cast<PlyProperty*>(&(_edge_props.at("vertex2"))));

    // setup face props
    describe_element_ply(ply_file, "face", _output_faces.size());
    describe_property_ply(ply_file, const_cast<PlyProperty*>(
                                        &(_face_props.at("vertex_indices"))));

    // add some comment
    append_comment_ply(ply_file, "Generated by VoxelCore");

    // header complete. write it.
    header_complete_ply(ply_file);

    // now start writing actual data
    put_element_setup_ply(ply_file, "vertex");
    for (int i = 0; i < _output_vts.size(); ++i)
    {
        put_element_ply(ply_file, const_cast<void*>(static_cast<const void*>(
                                      &_output_vts[i])));
    }
    put_element_setup_ply(ply_file, "edge");
    for (int i = 0; i < _output_edges.size(); ++i)
    {
        put_element_ply(ply_file, const_cast<void*>(static_cast<const void*>(
                                      &_output_edges[i])));
    }
    put_element_setup_ply(ply_file, "face");
    for (int i = 0; i < _output_faces.size(); ++i)
    {
        put_element_ply(ply_file, const_cast<void*>(static_cast<const void*>(
                                      &_output_faces[i])));
    }

    // that's it
    close_ply(ply_file);
    free_ply(ply_file);

    return ErrCode::SUCCESS;
}
} // namespace ply