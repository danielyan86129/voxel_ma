#include <voxelcore/spaceinfo.h>

const point SpaceConverter::vox_corner_offset[8] = {
    // x = -0.5 plane
    point(-0.5f, -0.5f, -0.5f), point(-0.5f, -0.5f, +0.5f),
    point(-0.5f, +0.5f, -0.5f), point(-0.5f, +0.5f, +0.5f),
    // x = +0.5 plane
    point(+0.5f, -0.5f, -0.5f), point(+0.5f, -0.5f, +0.5f),
    point(+0.5f, +0.5f, -0.5f), point(+0.5f, +0.5f, +0.5f)};