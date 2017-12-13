#ifndef COMMON_DEFS_H
#define COMMON_DEFS_H

#include <iostream>
#include <unordered_set>
#include <chrono>

#include <TriMesh.h>
#include <XForm.h>

using std::cout;
using std::endl;
using trimesh::point;
using trimesh::ivec3;
using trimesh::ivec4;
using trimesh::ivec2;
using trimesh::xform;
typedef trimesh::Vec<3, unsigned int> uTriFace;
//typedef trimesh::TriMesh::Face uTriFace;
typedef trimesh::vec3 TriColor;
struct RGBColor;

//#include <Eigen/Dense>
//using point = Eigen::Vector3f;
//using ivec3 = Eigen::Vector3i;
//using ivec4 = Eigen::Vector4i;
//using face = Eigen::Vector3i;

// the hasher for ivec2
struct ivec2Hash
{
	std::size_t operator()( ivec2 const& _t ) const
	{
		size_t h = 0;
		auto hasher = std::hash<int>();
		for ( size_t i = 0; i < 2; ++i )
		{
			// doing this: h = h*31 + hasher(_t[i])
			h = ( ( h << 5 ) - h ) + hasher( _t[ i ] );
		}
		return h; // or use boost::hash_combine
	}
};
// the hasher for ivec3
struct ivec3Hash
{
	std::size_t operator()( ivec3 const& _t ) const
	{
		size_t h = 0;
		auto hasher = std::hash<int>();
		for ( size_t i = 0; i < 3; ++i )
		{
			// doing this: h = h*31 + hasher(_t[i])
			h = ( ( h << 5 ) - h ) + hasher( _t[ i ] );
		}
		return h; // or use boost::hash_combine
	}
};

// simple struct for euler characteristic info
struct eulerchar
{
	int V, E, F, C, T, euler;
	eulerchar() {}
	eulerchar( int _v, int _e, int _f, int _t, int _c, int _euler )
	{
		V = _v;
		E = _e;
		F = _f;
		T = _t;
		C = _c;
		euler = _euler;
	}
	// print euler info to the std::cout
	inline void logToConsole( const char* _name_string )
	{
		cout << endl;
		cout << "------ Euler Characteristic (of " << _name_string << ") -------" << endl;
		cout << "euler = V - E + F = " << euler << endl;
		cout << "V / E / F / T / C = "
			<< V << " / " << E << " / "
			<< F << " / " << T << " / "
			<< C << endl;
		cout << "------ (Euler) -------" << endl;
	}
};

// rgb color structure
struct RGBColor
{
public:
	typedef unsigned char scalar_type;
	typedef trimesh::Vec<3, scalar_type> vec_type;
	RGBColor() 
	{
		m_c = vec_type( (scalar_type)255 );
	}
	explicit RGBColor( const TriColor& _c )
	{
		auto c = _c*255.0f;
		m_c[ 0 ] = std::max( std::min( c[ 0 ], 255.0f ), 0.0f );
		m_c[ 1 ] = std::max( std::min( c[ 1 ], 255.0f ), 0.0f );
		m_c[ 2 ] = std::max( std::min( c[ 2 ], 255.0f ), 0.0f );
	}
	scalar_type r() { return m_c[ 0 ]; }
	scalar_type g() { return m_c[ 1 ]; }
	scalar_type b() { return m_c[ 2 ]; }
private:
	vec_type m_c;
};

/// a light time-struct
typedef std::chrono::high_resolution_clock myclock;
struct timer
{
	// type representing any point of time
	typedef std::chrono::time_point<std::chrono::high_resolution_clock> timepoint;

	// type of duration in nano-sec
	typedef std::chrono::nanoseconds duration;
	typedef std::chrono::milliseconds millisec;

	// start and end timepoint
	timepoint s, e;
	// actual duration/time elapsed
	duration d;

	// ctor
	timer()
	{
		m_accum_time = false;
	}
	// start timer from a fresh state
	inline const timepoint& start()
	{
		s = myclock::now();
		d = duration::zero();
		return s;
	}
	// use this function if you want to accumulate chunks of time
	inline const timepoint& restart()
	{
		s = myclock::now();
		m_accum_time = true;
		return s;
	}
	inline const timepoint& stop()
	{
		e = myclock::now();
		if ( m_accum_time )
		{
			d += e - s;
			m_accum_time = false;
		}
		else
			elapse();
		return e;
	}
	inline duration elapse()
	{
		d = e - s;
		return d;
	}
	inline std::chrono::milliseconds elapseMilli() const
	{
		return std::chrono::duration_cast<millisec>( e - s );
	}

private:
	// when true, stop() will accumulate current duration 
	// (since last restart()) to total duration
	bool m_accum_time;
};

#endif // !COMMON_DEFS_H