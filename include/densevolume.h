#ifndef DENSEVOLUME_H
#define DENSEVOLUME_H

#include "Volume3DScalar.h"
#include <isosurface/volume.h> // Tao's volume rep.
#include <memory>

using std::shared_ptr;

//////////////////////////////////////////////////////////////////////////
/// This is the 3d dense volume representation
//////////////////////////////////////////////////////////////////////////
class DenseVolume : public Volume3DScalar
{
public:
    /*Constructors and destructor*/
    DenseVolume() {}
    DenseVolume(const shared_ptr<Volume>& _vol) : m_vol(_vol)
    {
        Volume3DScalar::m_sizes[0] = m_vol->getSizeX();
        Volume3DScalar::m_sizes[1] = m_vol->getSizeY();
        Volume3DScalar::m_sizes[2] = m_vol->getSizeZ();
    }
    ~DenseVolume() {}

    /*interfaces unique to DenseVolume*/
    // set the internal volume from a Tao's volume rep
    void setVol(const shared_ptr<Volume>& _vol)
    {
        m_vol = _vol;
        Volume3DScalar::m_sizes[0] = m_vol->getSizeX();
        Volume3DScalar::m_sizes[1] = m_vol->getSizeY();
        Volume3DScalar::m_sizes[2] = m_vol->getSizeZ();
    }
    /*implement parent interfaces*/
    double getDataAt(int _x, int _y, int _z)
    {
        return m_vol->getDataAt(_x, _y, _z);
    }
    double getDataAt(int _x, int _y, int _z) const
    {
        return m_vol->getDataAt(_x, _y, _z);
    }
    /*major interfaces consistent with Tao's definition*/
private:
    // pointer to Tao's volume instance
    shared_ptr<Volume> m_vol;
};

#endif