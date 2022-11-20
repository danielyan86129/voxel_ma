#pragma once
namespace voxelvoro
{

inline bool VoroInfo::isVertexValid(int _vi) const
{
    return onlyHasInside() || (m_vts_valid[_vi] && isVertexFinite(_vi));
}
inline bool VoroInfo::isEdgeValid(int _ei) const
{
    return onlyHasInside() || m_edge_valid[_ei];
}
inline bool VoroInfo::isFaceValid(int _fi) const
{
    return onlyHasInside() || m_face_valid[_fi];
}
inline bool VoroInfo::computeEdgeValidity(int _ei) const
{
    const auto e = m_geom.getEdge(_ei);
    return computeEdgeValidity(e);
}
inline bool VoroInfo::computeEdgeValidity(const ivec2& _e) const
{
    return isVertexValid(_e[0]) && isVertexValid(_e[1]);
}
inline bool VoroInfo::computeFaceValidity(int _fi) const
{
    vector<int> f_vts;
    m_geom.getFaceVRep(_fi, f_vts);
    for (auto v : f_vts)
        if (!isVertexValid(v))
            return false;
    return true;
}
inline bool VoroInfo::computeFaceValidity(const vector<int>& _f_vts) const
{
    for (auto v : _f_vts)
        if (!isVertexValid(v))
            return false;
    return true;
}
inline bool VoroInfo::computeFaceValidity(const uTriFace& _tri) const
{
    for (auto v : _tri)
        if (!isVertexValid(v))
            return false;
    return true;
}
inline bool VoroInfo::isVertexFinite(int _vi) const
{
    return onlyHasInside() || m_is_finite_v[_vi];
}
inline bool VoroInfo::isEdgeFinite(int _ei) const
{
    return onlyHasInside() || m_is_finite_e[_ei];
}
inline bool VoroInfo::isFaceFinite(int _fi) const
{
    return onlyHasInside() || m_is_finite_f[_fi];
}
inline bool VoroInfo::onlyHasInside() const { return m_only_inside > 0; }
inline bool VoroInfo::isMsureReady() const
{
    return m_v_msure.size() == geom().numVts() &&
           m_e_msure.size() == geom().numEdges() &&
           m_f_msure.size() == geom().numFaces();
}

inline ivec2 voxelvoro::VoroInfo::getSitesOfFace(int _fi) const
{
    return m_face_sites[_fi];
}
} // namespace voxelvoro