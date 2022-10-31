#include "octree.h"
#include "geomalgo.h"
#include "spaceinfo.h"
#include <cassert>
#include <filesystem>
#include <fstream>
#include <functional>
#include <unordered_set>

namespace fs = std::experimental::filesystem;
using std::ifstream;
using std::ofstream;
using std::unordered_set;

OctreeVolume::OctreeVolume() : Volume3DScalar()
{
    Volume3DScalar::m_sizes = ivec3(0);
    m_root = nullptr;
}

OctreeVolume::OctreeVolume(const char* _filename) : Volume3DScalar()
{
    m_data = nullptr;

    bool suc = false;
    // TODO: read entire file to buffer first, then parse the buffer.
    if (fs::path(_filename).extension() == ".sof")
        suc = read_sof_file(_filename);
    else if (fs::path(_filename).extension() == ".sog")
        suc = read_sog_file(_filename);
    if (suc)
        m_vol_filename = _filename;
    else
    {
        cout << "Octree volume construction failed." << endl;
    }
    // test
    /*auto test_ret = test_getDataAt();
    if ( test_ret == false )
    {
            cout << "Failed: test getDataAt() on all nodes." << endl;
    }
    else
    {
            cout << "Succeeded: test getDataAt() on all nodes." << endl;
    }*/
}

OctreeVolume::OctreeVolume(const Volume3DScalar* _src_vol)
{
    // max-node-res must be of 2^n
    m_max_res = std::max(std::max(_src_vol->getSizeX(), _src_vol->getSizeY()),
                         _src_vol->getSizeZ()) -
                1;
    m_max_res = 1 << util::log2Int(m_max_res);
    cout << "m_max_res = " << m_max_res << endl;
    // size of voxels (each sitting at corners of a node)
    Volume3DScalar::m_sizes =
        (_src_vol->getSizeX(), _src_vol->getSizeY(), _src_vol->getSizeZ());
    // actually build the tree
    m_root = OctreeVolume::build({0, 0, 0}, m_max_res, _src_vol);
}

OctreeVolume::~OctreeVolume()
{
    /*if ( m_data )
            delete[] m_data;*/
    if (m_root)
    {
        delnode(m_root);
        m_root = nullptr;
    }
}

bool OctreeVolume::writeToFile(const char* _file) const
{
    auto ext = fs::path(_file).extension().string();
    if (ext == ".sog")
        return write_sog_file(_file);
    else if (ext == ".sof")
        return write_sof_file(_file);
    else
    {
        printf("Octree file format (%s) not supported! \n", ext.c_str());
        return false;
    }
}

unsigned char OctreeVolume::getNodeValues(int _i, int _j, int _k) const
{
    unsigned char val;
    if (outOfVolume(_i, _j, _k))
    {
        val = 0;
    }
    else
    {
        auto data = m_data;
        auto retcode =
            get_node_values(_i, _j, _k, m_root, ivec3(0, 0, 0), m_max_res, val);
        assert(retcode == true);
    }
    return val;
}

bool OctreeVolume::outOfVolume(int _i, int _j, int _k) const
{
    return _i >= m_max_res || _j >= m_max_res || _k >= m_max_res;
}

void OctreeVolume::getBoundaryVoxels(vector<ivec3>& _voxels) const
{
    unordered_set<ivec3, ivec3Hash> voxels_set;
    int leaf_cnt = 0, internal_cnt = 0, empty_cnt = 0;
    int all_zero_leaf_cnt = 0, all_one_leaf_cnt = 0;

    // recurse until leaf node is reached,
    // then append all corners of that node to the list.
    std::function<void(OctreeNode*, ivec3, int)> grab_all_leaf =
        [&](OctreeNode* _node, ivec3 _off, int _len) {
            // get prepared if we are going to recurse
            unsigned char type;
            ivec3 node_off;
            int node_len = _len / 2;
            unsigned char values;

            // get type of the node
            type = _node->type;

            // only care about internal and leaf nodes
            if (type == 2)
            {
                leaf_cnt++;
                // leaf node. no need to recurse.
                auto leafnode = dynamic_cast<LeafNode*>(_node);
                values = leafnode->getValues();
                if (values == 0)
                {
                    all_zero_leaf_cnt++;
                }
                else if (values == 0xFF)
                {
                    all_one_leaf_cnt++;
                }
                // add all its corners as boundary voxels
                for (int i = 0; i < 2; i++)
                    for (int j = 0; j < 2; j++)
                        for (int k = 0; k < 2; k++)
                        {
                            voxels_set.insert(_off + ivec3(k, j, i));
                        }
            }
            else if (type == 0)
            {
                internal_cnt++;
                // interior node. recurse.
                auto intnode = dynamic_cast<IntNode*>(_node);
                for (int i = 0; i < 2; i++)
                    for (int j = 0; j < 2; j++)
                        for (int k = 0; k < 2; k++)
                        {
                            node_off =
                                _off + ivec3(/*k, j, i*/ i, j, k) * node_len;
                            auto childnode =
                                intnode->getChild(/*k, j, i*/ i, j, k);
                            grab_all_leaf(childnode, node_off, node_len);
                        }
            }
        }; // grab_all_leaf()

    grab_all_leaf(m_root, ivec3(0, 0, 0), m_max_res);
    cout << "# leaves with all 0s/1s: " << all_zero_leaf_cnt << "/"
         << all_one_leaf_cnt << endl;
    /*cout << "# leaf/empty/internal nodes: "
            << leaf_cnt << "/" << empty_cnt << "/" << internal_cnt << endl;*/

    // copy all unique boundary voxels to output list
    _voxels.clear();
    _voxels.reserve(voxels_set.size());
    for (auto it = voxels_set.begin(); it != voxels_set.end();)
    {
        _voxels.push_back(*it);
        it = voxels_set.erase(it);
    }
}

trimesh::Box<3, int> OctreeVolume::getBoundingBox() const
{
    trimesh::Box<3, int> bb;
    vector<ivec3> bndry_voxels;
    getBoundaryVoxels(bndry_voxels);
    for (const auto& v : bndry_voxels)
    {
        bb += v;
    }
    return bb;
}

double OctreeVolume::getDataAt(int _x, int _y, int _z)
{
    // cout << "finding voxel: " << ivec3( _x, _y, _z ) << " ";
    //  if any of x, y, or z is out-of-range from below,
    //  simply return 0
    /*if ( _x < 0 || _y < 0 || _z < 0 )
            return 0;*/
    // else, one or more of them may be out-of-range from above
    // in this case, we only need to clamp it back into the range
    // if it is equal to the higher end of the range
    unsigned char offset = 0;
    if (_x == m_max_res)
    {
        _x = m_max_res - 1;
        offset |= 4; // set 00000100
    }
    if (_y == m_max_res)
    {
        _y = m_max_res - 1;
        offset |= 2; // set 00000010
    }
    if (_z == m_max_res)
    {
        _z = m_max_res - 1;
        offset |= 1; // set 00000001
    }
    auto val = getNodeValues(_x, _y, _z);
    // cout << endl;
    //  return the the lowest left voxel
    auto ret = (((val >> offset) & 1) == 0) ? 1 : -1;
    return ret;
}

double OctreeVolume::getDataAt(int _x, int _y, int _z) const
{
    unsigned char offset = 0;
    if (_x == m_max_res)
    {
        _x = m_max_res - 1;
        offset |= 4; // 00000100
    }
    if (_y == m_max_res)
    {
        _y = m_max_res - 1;
        offset |= 2; // 00000010
    }
    if (_z == m_max_res)
    {
        _z = m_max_res - 1;
        offset |= 1; // 00000001
    }
    auto val = getNodeValues(_x, _y, _z);
    // return the the lowest left voxel
    auto ret = (((val >> offset) & 1) == 0) ? 1 : -1;
    return ret;
}

// construct octree as data being read from .sof file
// .sof has no transformation info, so transform mat stays identity
bool OctreeVolume::read_sof_file(const string& _sof_file)
{
    ifstream infile(_sof_file, std::ios::binary);
    if (!infile.is_open())
    {
        cout << "Error: cannot open .sof file: " << _sof_file << endl;
        return false;
    }

    // process header & init size info
    infile.read((char*)(&m_max_res), sizeof(int));
    Volume3DScalar::m_sizes = ivec3(m_max_res + 1);
    cout << "octree res = " << m_max_res << endl;
    Volume3DScalar::m_vox_to_mesh = trimesh::xform::identity();
    Volume3DScalar::m_mesh_to_vox = trimesh::xform::identity();

    // recursively create octree from the rest of the file
    int leaf_cnt = 0, internal_cnt = 0, empty_cnt = 0;
    std::function<OctreeNode*(ivec3, int)> mknode =
        [&](const ivec3& _off, int _len) -> OctreeNode* {
        // get prepared if we are going to recurse
        unsigned char type;
        ivec3 node_off;
        int node_len = _len / 2;
        unsigned char values;

        // get type of the node
        infile.read((char*)&type, 1);

        if (type == LEAF)
        {
            leaf_cnt++;
            // leaf node. no need to recurse.
            infile.read((char*)&values, 1);
            auto nd = new LeafNode(values);

            return nd;
        }
        else if (type == INTERNAL)
        {
            internal_cnt++;
            // interior node. recurse.
            auto nd = new IntNode();
            OctreeNode* cnode;
            // note: in the file, the children of internal node
            // are stored in this order
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    for (int k = 0; k < 2; k++)
                    {
                        node_off = _off + node_len * ivec3(i, j, k);
                        cnode = mknode(node_off, node_len);
                        nd->setChild(i, j, k, cnode);
                    }

            return nd;
        }
        else if (type == EMPTY)
        {
            empty_cnt++;
            // empty node. skip its sign value.
            infile.read((char*)&values, 1);
            auto nd = new EmptyNode(values);

            /*if ( values == 0 )
            cout << "inside empty node found. " << endl;*/

            return nd;
        }
        else
        {
            // This should never happen!
            cout << "Error: invalid node type " << (int)type << "!" << endl;
            cout << "octree so far constructed - Leaf/empty/internal: "
                 << leaf_cnt << "/" << empty_cnt << "/" << internal_cnt << endl
                 << endl;
            exit(-1);
        }
    };

    cout << "------- creating octree -------- " << endl;
    m_root = mknode(ivec3(0, 0, 0), m_max_res);
    cout << "octree constructed. Leaf/empty/internal: " << leaf_cnt << "/"
         << empty_cnt << "/" << internal_cnt << endl
         << endl;

    // done.
    infile.close();
    return true;
}

// construct octree as data being read from .sog file
// transform mat will reflect the transformation info in the file
bool OctreeVolume::read_sog_file(const string& _sog_file)
{
    ifstream infile(_sog_file, std::ios::binary);
    if (!infile.is_open())
    {
        cout << "Error: cannot open .sog file: " << _sog_file << endl;
        return false;
    }

    // process header & init size info
    auto titlestr = "SOG.Format 1.0";
    auto header_buf = new char[128];
    auto buf = header_buf;
    infile.read(buf, 128);
    buf += std::strlen(titlestr) + 1; // skip title + null terminator
    float x, y, z;
    x = *(float*)(buf);
    buf += sizeof(float);
    y = *(float*)(buf);
    buf += sizeof(float);
    z = *(float*)(buf);
    buf += sizeof(float);
    auto lowerleft = point(x, y, z);
    float bbox_len_input_space = *(float*)(buf);
    buf += sizeof(float);
    delete[] header_buf;

    Volume3DScalar::m_vox_to_mesh = trimesh::xform::identity();
    Volume3DScalar::m_mesh_to_vox = trimesh::xform::identity();

    // recursively create octree from the rest of the file
    int leaf_cnt = 0, internal_cnt = 0, empty_cnt = 0;
    int n_levels = 0;
    float dump[3];
    // start from root at level 0, recursively construct the tree
    std::function<OctreeNode*(int /*ivec3, int*/)> mknode =
        [&](int _level /*const ivec3& _off, int _len*/) -> OctreeNode* {
        // get prepared if we are going to recurse
        unsigned char type;
        /*ivec3 node_off;
        int node_len = _len / 2;*/
        unsigned char values;

        // get type of the node
        infile.read((char*)&type, 1);

        if (type == LEAF)
        {
            n_levels = std::max(n_levels, _level);
            leaf_cnt++;
            // leaf node. no need to recurse.
            infile.read((char*)&values, 1);
            auto nd = new LeafNode(values);
            // depelete rest geometry info (3 floats)
            infile.read((char*)&dump, 3 * sizeof(float));

            return nd;
        }
        else if (type == INTERNAL)
        {
            internal_cnt++;
            // interior node. recurse.
            auto nd = new IntNode();
            OctreeNode* cnode;
            // note: in the file, the children of internal node
            // are stored in this order
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    for (int k = 0; k < 2; k++)
                    {
                        /*node_off = _off + node_len * ivec3( i, j, k );*/
                        cnode = mknode(_level + 1 /*node_off, node_len*/);
                        nd->setChild(i, j, k, cnode);
                    }

            return nd;
        }
        else if (type == EMPTY)
        {
            empty_cnt++;
            // empty node. skip its sign value.
            infile.read((char*)&values, 1);
            auto nd = new EmptyNode(values);

            /*if ( values == 0 )
            cout << "inside empty node found. " << endl;*/

            return nd;
        }
        else
        {
            // This should never happen!
            cout << "Error: invalid node type " << (int)type << "!" << endl;
            cout << "octree so far constructed - Leaf/empty/internal: "
                 << leaf_cnt << "/" << empty_cnt << "/" << internal_cnt << endl
                 << endl;
            exit(-1);
        }
    };

    cout << "------- creating octree -------- " << endl;
    infile.read((char*)&m_max_res, sizeof(int));
    cout << "octree res = " << m_max_res << endl;
    m_root = mknode(0 /*ivec3( 0, 0, 0 ), m_max_res*/);
    // m_max_res = std::pow( 2, n_levels );
    // cout << "octree res = " << m_max_res << endl;
    Volume3DScalar::m_sizes = ivec3(m_max_res + 1);

    auto s = (double)bbox_len_input_space / (double)m_max_res;
    auto scale_mat = trimesh::xform::scale(s);
    auto trans_mat = trimesh::xform::trans(lowerleft);
    Volume3DScalar::m_vox_to_mesh = trans_mat * scale_mat;
    cout << "bbox len in mesh space: " << bbox_len_input_space << endl;
    cout << "scale factor: " << s << endl;
    cout << "lower left: " << lowerleft << endl;
    // cout << "vox to mesh space transformation: " <<
    // Volume3DScalar::m_vox_to_mesh << endl;
    //  inverse transformation
    scale_mat = trimesh::xform::scale(1.0 / s);
    trans_mat = trimesh::xform::trans(-lowerleft);
    Volume3DScalar::m_mesh_to_vox = scale_mat * trans_mat;
    auto auto_mesh_to_vox = trimesh::inv(Volume3DScalar::m_vox_to_mesh);
    for (auto i = 0; i < auto_mesh_to_vox.size(); ++i)
    {
        auto v1 = auto_mesh_to_vox[i];
        auto v2 = Volume3DScalar::m_vox_to_mesh[i];
        if (!util::is_equal(v1, v2, bbox_len_input_space * 0.0000001))
        {
            cout
                << "Warning: mesh_to_vox and auto_vox_to_mesh are not the same!"
                << endl;
            cout << "mesh-to-vox: " << Volume3DScalar::m_mesh_to_vox << endl;
            cout << "auto-mesh-to-vox: " << auto_mesh_to_vox << endl;
            break;
        }
    }
    cout << "octree constructed. Leaf/empty/internal: " << leaf_cnt << "/"
         << empty_cnt << "/" << internal_cnt << endl
         << endl;

    // done.
    infile.close();
    return true;
}

bool OctreeVolume::write_sof_file(const string& _sof_file) const
{
    int n_leaf = 0, n_empty = 0, n_internal = 0;
    ofstream os(_sof_file, std::ios::binary);
    if (!os.is_open())
    {
        cout << "Error: cannot write octree to file " << _sof_file << endl;
        return false;
    }
    // write node res.
    cout << "writing max-res = " << m_max_res << endl;
    os.write((char*)&m_max_res, sizeof(int));

    std::function<void(const OctreeNode*)> write_node =
        [&](const OctreeNode* _node) -> bool {
        auto type = (char)_node->getType();
        os.write(&type, sizeof(char));
        unsigned char values = 0x00;
        if (type == NodeType::LEAF)
        {
            n_leaf++;
            auto nd = dynamic_cast<const LeafNode*>(_node);
            values = nd->getValues();
            os.write((char*)&values, sizeof(values));
            /*if ( values == 0xa )
            {
            int stop = 1;
            os.close();
            exit( -2 );
            }*/
        }
        else if (type == NodeType::EMPTY)
        {
            n_empty++;
            auto nd = dynamic_cast<const EmptyNode*>(_node);
            values = nd->getValues();
            os.write((char*)&values, sizeof(values));
            /*if ( values == 0xa )
            {
            int stop = 1;
            os.close();
            exit( -2 );
            }*/
        }
        else if (type == NodeType::INTERNAL)
        {
            n_internal++;
            auto nd = dynamic_cast<const IntNode*>(_node);
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    for (int k = 0; k < 2; k++)
                        write_node(nd->getChild(i, j, k));
        }
        else
        {
            // This should never happen!
            cout << "Error: invalid node type " << (int)type << "!" << endl;
            return false;
        }
        return true;
    };

    // write everything hanging from root
    write_node(m_root);
    os.close();
    cout << "# Octree nodes (leaf/empty/internal): " << n_leaf << "/" << n_empty
         << "/" << n_internal << endl;
    return true;
}
bool OctreeVolume::write_sog_file(const string& _sog_file) const
{
    int n_leaf = 0, n_empty = 0, n_internal = 0;
    ofstream os(_sog_file, std::ios::binary);
    if (!os.is_open())
    {
        cout << "Error: cannot write octree to file " << _sog_file << endl;
        return false;
    }

    /* prepare header */
    int header_size = 128;
    char* header_buf = new char[header_size];
    auto buf = header_buf;
    // starts with title
    const char* title = "SOG.Format 1.0";
    int title_len = std::strlen(title) + 1;
    memcpy(buf, title, title_len);
    buf += title_len;
    // lower-left. TODO: find a way to extract this info from vox-mesh-matrix
    point vol_o(0.0f);
    int elem_size = sizeof(float);
    memcpy(buf, (char*)&vol_o[0], elem_size);
    buf += elem_size;
    memcpy(buf, (char*)&vol_o[0], elem_size);
    buf += elem_size;
    memcpy(buf, (char*)&vol_o[0], elem_size);
    buf += elem_size;
    // vol len
    float vol_l = this->m_max_res;
    memcpy(buf, (char*)&vol_l, sizeof(vol_l));
    buf += sizeof(float);

    /*write header*/
    cout << "writing header (with max-res: " << m_max_res << ")" << endl;
    os.write(header_buf, header_size);
    // finally write max-res
    os.write((char*)&this->m_max_res, sizeof(int));

    /* cleanup space */
    delete header_buf;

    std::function<void(const OctreeNode*, const ivec3&, int)> write_node =
        [&](const OctreeNode* _node, const ivec3& _off, int _len) -> bool {
        auto type = (char)_node->getType();
        os.write(&type, sizeof(char));
        unsigned char values = 0x00;
        if (type == NodeType::LEAF)
        {
            n_leaf++;
            auto nd = dynamic_cast<const LeafNode*>(_node);
            values = nd->getValues();
            os.write((char*)&values, sizeof(values));
            // TODO: now use center as feature. consider use better geometry.
            point feat_pt(_off[0] + 0.5, _off[1] + 0.5, _off[2] + 0.5);
            os.write((char*)(feat_pt.data()), sizeof(float) * 3);
            /*if ( values == 0xa )
            {
            int stop = 1;
            os.close();
            exit( -2 );
            }*/
        }
        else if (type == NodeType::EMPTY)
        {
            n_empty++;
            auto nd = dynamic_cast<const EmptyNode*>(_node);
            values = nd->getValues();
            os.write((char*)&values, sizeof(values));
            /*if ( values == 0xa )
            {
            int stop = 1;
            os.close();
            exit( -2 );
            }*/
        }
        else if (type == NodeType::INTERNAL)
        {
            n_internal++;
            auto nd = dynamic_cast<const IntNode*>(_node);
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    for (int k = 0; k < 2; k++)
                    {
                        int c_len = _len / 2;
                        ivec3 c_off = _off + ivec3(i, j, k) * c_len;
                        write_node(nd->getChild(i, j, k), c_off, c_len);
                    }
        }
        else
        {
            // This should never happen!
            cout << "Error: invalid node type " << (int)type << "!" << endl;
            return false;
        }
        return true;
    }; // write_node()

    // write everything hanging from root
    write_node(m_root, ivec3(0), m_max_res);
    os.close();
    cout << "# Octree nodes (leaf/empty/internal): " << n_leaf << "/" << n_empty
         << "/" << n_internal << endl;
    return true;
}

bool OctreeVolume::get_node_values(int _i, int _j, int _k, OctreeNode* _node,
                                   ivec3 _off, int _len,
                                   unsigned char& _vals) const
{
    bool found = false;
    unsigned char ret_values;

    // get prepared if we are going to recurse
    char type;
    ivec3 node_off;
    int node_len = _len / 2;

    // get type of the node
    type = _node->type;

    // see if the requested node is within the region
    // return immediately if not
    ivec3 ncoord(_i, _j, _k);
    auto ncoord_t = ncoord - _off;
    if (ncoord_t[0] < 0 || ncoord_t[0] >= _len || ncoord_t[1] < 0 ||
        ncoord_t[1] >= _len || ncoord_t[2] < 0 || ncoord_t[2] >= _len)
    {
        return false;
    }

    // cout << "-> ";
    // if ( type == 2 )
    //	cout << "lf ";
    // else if ( type == 1 )
    //	cout << "et ";
    // else
    //	cout << "in ";
    // cout << "{" << _off << "," << _len << "} "; // debug

    if (type == 2)
    {
        if (_off != ivec3(_i, _j, _k))
        {
            cout << "Assertion off == (i,j,k) failed! " << _off << ", "
                 << ivec3(_i, _j, _k) << endl;
            exit(-1);
        }
        // leaf node. no need to recurse.
        auto leafnode = dynamic_cast<LeafNode*>(_node);
        found = true;
        ret_values = leafnode->getValues();
        // cout << "leaf node value: " << (int)ret_values << endl;
        /*if ( ret_values == 0 )
                cout << "(inside) leaf node value: " << (int)ret_values << endl;
        else if ( ret_values > 0 )
                cout << "(boundary) leaf node value: " << (int)ret_values <<
        endl;*/
    }
    else if (type == 1)
    {
        // empty node. no need to recurse.
        auto emptynode = dynamic_cast<EmptyNode*>(_node);
        found = true;
        // replicate this value to all 8 bits and return
        if (emptynode->getValues() == 1)
            ret_values = 0xFF;
        else
            ret_values = 0;

        /*if ( ret_values == 0 )
                cout << "empty node value: " << (int)ret_values << endl;*/
    }
    else if (type == 0)
    {
        // interior node. recurse.
        auto intnode = dynamic_cast<IntNode*>(_node);
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                for (int k = 0; k < 2; k++)
                {
                    node_off = _off + ivec3(/*k, j, i*/ i, j, k) * node_len;
                    unsigned char node_values;
                    auto childnode = intnode->getChild(/*k, j, i*/ i, j, k);
                    found = get_node_values(_i, _j, _k, childnode, node_off,
                                            node_len, node_values);
                    if (found)
                    {
                        ret_values = node_values;
                        goto RETURN;
                    }
                }
    }
    else
    {
        // This should never happen!
        cout << "Error: invalid node type " << type << "!" << endl;
        found = false;
    }

RETURN:
    _vals = ret_values;
    return found;
}

bool OctreeVolume::test_getDataAt()
{
    bool ret_code = false;
    ret_code = test_node(m_root, ivec3(0, 0, 0), m_max_res);
    return ret_code;
}

bool OctreeVolume::test_node(OctreeNode* _node, ivec3 _off, int _len)
{
    bool ret_code = true;
    auto type = _node->type;
    ivec3 vox;
    int value1, value2;

    if (type == 2)
    {
        auto leafnode = dynamic_cast<LeafNode*>(_node);
        auto nodevalue = leafnode->getValues();
        // test 8 corners with the value returned by getDataAt()
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                for (int k = 0; k < 2; ++k)
                {
                    // the value of corner returned by getDataAt()
                    vox = ivec3(/*k, j, i*/ i, j, k) + _off;
                    value1 = (int)getDataAt(vox[0], vox[1], vox[2]);
                    // get the value of corner from this node directly
                    value2 = leafnode->getCorner(/*k, j, i*/ i, j, k);
                    value2 = value2 == 0 ? 1 : -1;

                    if (value1 != value2)
                    {
                        cout << "values doesn't match! At node {" << (int)type
                             << ", " << _off << ", " << _len << "}, iter "
                             << ivec3(i, j, k) << ". value1/value2 = " << value1
                             << "/" << value2 << endl;
                        ret_code = false;
                    }
                }
    }
    else if (type == 1)
    {
        // empty node
        auto emptynode = dynamic_cast<EmptyNode*>(_node);
        auto nodevalue = emptynode->getValues();
        for (int i = 0; i <= _len; ++i)
            for (int j = 0; j <= _len; ++j)
                for (int k = 0; k <= _len; ++k)
                {
                    vox = _off + ivec3(/*k, j, i*/ i, j, k);
                    value1 = (int)getDataAt(vox[0], vox[1], vox[2]);
                    value2 = nodevalue;
                    value2 = value2 == 0 ? 1 : -1;

                    if (value1 != value2)
                    {
                        cout << "values doesn't match! {" << (int)type << ", "
                             << _off << ", " << _len << "}, iter "
                             << ivec3(i, j, k) << ". value1/value2 = " << value1
                             << "/" << value2 << endl;
                        ret_code = false;
                    }
                }
    }
    else
    {
        auto intnode = dynamic_cast<IntNode*>(_node);
        ivec3 noff;
        int nlen = _len / 2;
        ;
        // recursively test 8 children
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                for (int k = 0; k < 2; ++k)
                {
                    noff = _off + ivec3(/*k, j, i*/ i, j, k) * nlen;
                    auto node_ret = test_node(
                        intnode->getChild(/*k, j, i*/ i, j, k), noff, nlen);
                    if (node_ret == false)
                    {
                        ret_code = false;
                    }
                }
    }

RETURN:
    return ret_code;
}

OctreeVolume::OctreeNode* OctreeVolume::build(const ivec3& _p, int _l,
                                              const Volume3DScalar* _src_vol)
{
    /* region (p, l) is outside of _src_vol completely*/
    if (!util::intersect(
            _p, {_l, _l, _l}, {0, 0, 0},
            {_src_vol->getSizeX(), _src_vol->getSizeY(), _src_vol->getSizeZ()}))
    {
        // cout << "empty node. out-of-bound." << endl;
        return new EmptyNode(0x01); // 1: outside, 0: inside
    }
    else
    {
        // cout << "intersected!" << endl;
        if (_l == 1) // bottom octree node. collect values of leaves below it.
        {
            auto node = new OctreeVolume::LeafNode(0x00);
            for (auto i = 0; i < 2; ++i)
                for (auto j = 0; j < 2; ++j)
                    for (auto k = 0; k < 2; ++k)
                    {
                        auto off = ivec3(/*k, j, i*/ i, j, k);
                        auto vox_val = SpaceConverter::get_occupancy_at_vox(
                            _p + off, _src_vol);
                        vox_val = vox_val == 1 ? 0 : 1;
                        // if (vox_val != 0) cout << "leaf corner value: " <<
                        // vox_val << endl;
                        node->setCorner(off[0], off[1], off[2], vox_val);
                    }
            // cout << "leaf node. values: " << (int)node->getValues() << endl;
            return node;
        }
        else // ( _l >= 2 ), node on interior levels. subdivide.
        {
            // cout << "internal node." << endl;
            auto node = new OctreeVolume::IntNode();
            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j)
                    for (int k = 0; k < 2; ++k)
                    {
                        auto childlen = _l / 2;
                        auto off = ivec3(/*k, j, i*/ i, j, k);
                        const auto childpos = _p + off * childlen;
                        auto child = build(childpos, childlen, _src_vol);
                        node->setChild(off[0], off[1], off[2], child);
                    }

            // check if need to merge children
            unsigned char val, cur_val;
            auto is_homo = node->getChild(0)->isHomo(val); // init correctly
            for (auto i = 0; i < 8; ++i)
            {
                is_homo = is_homo && node->getChild(i)->isHomo(cur_val);
                if (!is_homo || cur_val != val)
                {
                    // cout << "not-merged. l=" << _l << endl;
                    return node;
                }
            }
            // cout << "merged. l=" << _l << endl;
            //  all children having same value.
            //  dealloc old node and recreate an EmptyNode to return
            OctreeVolume::delnode(node);
            return new OctreeVolume::EmptyNode(val);
        }
    }
}

void OctreeVolume::walkTree(Walker* _w)
{
    walkSubtree(_w, m_root, ivec3(0, 0, 0), m_max_res);
}

void OctreeVolume::walkSubtree(Walker* _w, OctreeNode* _node, const ivec3& _off,
                               int _len)
{
    if (_node->getType() == LEAF || _node->getType() == EMPTY)
    {
        if (_node->getType() == LEAF)
        {
            _w->leaf(this, _node, _off, _len);
        }
        return;
    }
    // if ( _node->getType() == EMPTY )
    //{
    //	//_w->empty();
    // }
    if (!_w->pre(this, _node, _off, _len))
        return;
    // recurse into this internal node
    auto intnode = dynamic_cast<IntNode*>(_node);
    ivec3 noff;
    int nlen;
    for (int ix = 0; ix < 2; ++ix)
        for (int iy = 0; iy < 2; ++iy)
            for (int iz = 0; iz < 2; ++iz)
            {
                auto cnode = intnode->getChild(ix, iy, iz);
                nlen = _len / 2;
                noff = _off + ivec3(ix, iy, iz) * nlen;
                walkSubtree(_w, cnode, noff, nlen);
            }
    _w->post(this, _node, _off, _len);
}
