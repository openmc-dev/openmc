/*
 * Copyright (c) 2005 Lawrence Livermore National Laboratory under
 * contract number B545069 with the University of Wisconsin - Madison.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <math.h>

#ifndef PYNE_IS_AMALGAMATED
#include "measure.h"
#endif

class CartVect {
  private:
    double coords[3];
    
  public:
  
    inline CartVect() {}
    
    inline CartVect( double tx, double ty, double tz ) { set(tx,ty,tz); }
    
    inline CartVect( const CartVect& other ) { set( other.coords); }
    
    inline void set( double tx, double ty, double tz )
      { coords[0] = tx; coords[1] = ty; coords[2] = tz; }
    
    inline void set( const double* c )
      { coords[0] = c[0]; coords[1] = c[1]; coords[2] = c[2]; }
    
    inline double x() const { return coords[0]; }
    inline double y() const { return coords[1]; }
    inline double z() const { return coords[2]; }
  
    inline CartVect& operator+=( const CartVect& other )
    {
      coords[0] += other.coords[0];
      coords[1] += other.coords[1];
      coords[2] += other.coords[2];
      return *this;
    }
    
    inline CartVect& operator-=( const CartVect& other )
    {
      coords[0] -= other.coords[0];
      coords[1] -= other.coords[1];
      coords[2] -= other.coords[2];
      return *this;
    }
    
    inline CartVect& operator*=( const CartVect& other );

    inline double lensqr() const;
    
    inline double len() const;
    
    inline CartVect operator~( ) const;

      
    inline CartVect& operator*=( double a )
    {
      coords[0] *= a;
      coords[1] *= a;
      coords[2] *= a;
      return *this;
    }
    
    inline CartVect& operator/=( double a )
    {
      coords[0] /= a;
      coords[1] /= a;
      coords[2] /= a;
      return *this;
    }
   
    
};

inline CartVect operator+( const CartVect& v1, const CartVect& v2 )
{
  CartVect rval(v1);
  rval += v2;
  return rval;
}

inline CartVect operator-( const CartVect& v1, const CartVect& v2 )
{
  CartVect rval(v1);
  rval -= v2;
  return rval;
}

inline double operator%( const CartVect& v1, const CartVect& v2 )
{
  return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
}

inline CartVect operator*( const CartVect& v1, const CartVect& v2 )
{
  return CartVect( v1.y() * v2.z() - v1.z() * v2.y(),
                   v1.z() * v2.x() - v1.x() * v2.z(),
                   v1.x() * v2.y() - v1.y() * v2.x() );
}

inline CartVect CartVect::operator~() const
{
  double invlen = 1.0 / len();
  return CartVect( invlen * x(), invlen * y(), invlen * z() );
}
     
inline CartVect& CartVect::operator*=( const CartVect& other )
      { return *this = *this * other; }

inline double CartVect::lensqr() const
      { return *this % *this; }
    
inline double CartVect::len() const
      { return sqrt(lensqr()); }
 
inline static double tet_volume( const CartVect& v0,
                                 const CartVect& v1,
                                 const CartVect& v2, 
                                 const CartVect& v3 )
{
  return 1./6. * ( ((v1 - v0) * (v2 - v0)) % (v3 - v0) );
}

double edge_length( const double* start_vtx_coords,
                    const double*   end_vtx_coords )
{
  const CartVect* start = reinterpret_cast<const CartVect*>(start_vtx_coords);
  const CartVect*   end = reinterpret_cast<const CartVect*>(  end_vtx_coords);
  return (*start - *end).len();
}

double measure( moab::EntityType type,
                int num_vertices,
                const double* vertex_coordinates )
{
  const CartVect* coords = reinterpret_cast<const CartVect*>(vertex_coordinates);
  switch( type )
  {
    case moab::MBEDGE:
      return (coords[0] - coords[1]).len();
    case moab::MBTRI:
      return 0.5 * ((coords[1] - coords[0]) * (coords[2] - coords[0])).len();
    case moab::MBQUAD:
      num_vertices = 4;
    case moab::MBPOLYGON:
    {
      CartVect mid(0,0,0);
      for (int i = 0; i < num_vertices; ++i)
        mid += coords[i];
      mid /= num_vertices;
      
      double sum = 0.0;
      for (int i = 0; i < num_vertices; ++i)
      {
        int j = (i+1)%num_vertices;
        sum += ((mid - coords[i]) * (mid - coords[j])).len();
      }
      return 0.5 * sum;
    }
    case moab::MBTET:
      return tet_volume( coords[0], coords[1], coords[2], coords[3] ) ;
    case moab::MBPYRAMID:
      return tet_volume( coords[0], coords[1], coords[2], coords[4] ) +
             tet_volume( coords[0], coords[2], coords[3], coords[4] ) ;
    case moab::MBPRISM:
      return tet_volume( coords[0], coords[1], coords[2], coords[5] ) +
             tet_volume( coords[3], coords[5], coords[4], coords[0] ) +
             tet_volume( coords[1], coords[4], coords[5], coords[0] ) ;
    case moab::MBHEX:
      return tet_volume( coords[0], coords[1], coords[3], coords[4] ) +
             tet_volume( coords[7], coords[3], coords[6], coords[4] ) +
             tet_volume( coords[4], coords[5], coords[1], coords[6] ) +
             tet_volume( coords[1], coords[6], coords[3], coords[4] ) +
             tet_volume( coords[2], coords[6], coords[3], coords[1] ) ;
    default:
      return 0.0;
  }
}
      
