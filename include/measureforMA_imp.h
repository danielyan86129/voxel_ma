inline float MeasureForMA::lambdaForFace(const point& _a, const point& _b)
{
    return trimesh::dist(_a, _b);
}