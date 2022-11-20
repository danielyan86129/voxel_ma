#pragma once

#include "commondefs.h"

/************************************************************************
 * compute various measure on medial axis
 ************************************************************************/

class MeasureForMA
{
public:
    enum meassuretype
    {
        LAMBDA
    };
    //
    // compute the lambda measure value for a face of a ma, given 2 closest
    // boundary points.
    //
    static inline float lambdaForFace(const point& _a, const point& _b);

private:
    MeasureForMA() = default;
    MeasureForMA(const MeasureForMA&) = default;
    MeasureForMA& operator=(const MeasureForMA&) = default;
};

#include "measureforMA_imp.h"