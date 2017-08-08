#include <cstring>  // For strcmp
#include <iostream>
#include <limits>  // For numeric_limits
#include <math.h>  // For fabs
#include "pugixml/pugixml.hpp"

//==============================================================================
// Constants
//==============================================================================

const double FP_COINCIDENT = 1e-12;
const double INFTY = std::numeric_limits<double>::max();

//==============================================================================
// Global array of surfaces
//==============================================================================

class Surface;
Surface **surfaces_c;

//==============================================================================
// SURFACE type defines a first- or second-order surface that can be used to
// construct closed volumes (cells)
//==============================================================================

class Surface
{
  int id;                    // Unique ID
  int neighbor_pos[],        // List of cells on positive side
      neighbor_neg[];        // List of cells on negative side
  int bc;                    // Boundary condition
  //TODO: swith that zero to a NONE constant.
  int i_periodic = 0;        // Index of corresponding periodic surface
  char name[104];            // User-defined name

public:
  bool sense(double xyz[3], double uvw[3]);
  void reflect(double xyz[3], double uvw[3]);
  virtual double evaluate(double xyz[3]) = 0;
  virtual double distance(double xyz[3], double uvw[3], bool coincident) = 0;
  virtual void normal(double xyz[3], double uvw[3]) = 0;
};

bool
Surface::sense(double xyz[3], double uvw[3])
{
  // Evaluate the surface equation at the particle's coordinates to determine
  // which side the particle is on.
  const double f = evaluate(xyz);

  // Check which side of surface the point is on.
  bool s;
  if (fabs(f) < FP_COINCIDENT)
  {
    // Particle may be coincident with this surface. To determine the sense, we
    // look at the direction of the particle relative to the surface normal (by
    // default in the positive direction) via their dot product.
    double norm[3];
    normal(xyz, norm);
    s = (uvw[0] * norm[0] + uvw[1] * norm[1] + uvw[2] * norm[2] > 0.0);
  }
  else
  {
    s = (f > 0.0);
  }
  return s;
}

void
Surface::reflect(double xyz[3], double uvw[3])
{

  // Determine projection of direction onto normal and squared magnitude of
  // normal.
  double norm[3];
  normal(xyz, norm);
  const double projection = norm[0]*uvw[0] + norm[1]*uvw[1] + norm[2]*uvw[2];
  const double magnitude = norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2];

  // Reflect direction according to normal.
  uvw[0] -= 2.0 * projection / magnitude * norm[0];
  uvw[1] -= 2.0 * projection / magnitude * norm[1];
  uvw[2] -= 2.0 * projection / magnitude * norm[2];
}

//==============================================================================

class SurfaceZPlane : public Surface
{
  double z0;
public:
  SurfaceZPlane(pugi::xml_node surf_node);
  double evaluate(double xyz[3]);
  double distance(double xyz[3], double uvw[3], bool coincident);
  void normal(double xyz[3], double uvw[3]);
};

SurfaceZPlane::SurfaceZPlane(pugi::xml_node surf_node)
{
  const char *coeffs = surf_node.attribute("coeffs").value();
  int stat = sscanf(coeffs, "%lf", &z0);
  if (stat != 1)
  {
    std::cout << "Something went wrong reading surface coeffs!" << std::endl;
  }
}

double
SurfaceZPlane::evaluate(double xyz[3])
{
  double f = xyz[2] - z0;
  return f;
}

double
SurfaceZPlane::distance(double xyz[3], double uvw[3], bool coincident)
{
  double f = z0 - xyz[2];
  double d;
  if (coincident or fabs(f) < FP_COINCIDENT or uvw[2] == 0.0)
  {
    d = INFTY;
  }
  else
  {
    d = f / uvw[2];
    if (d < 0.0) d = INFTY;
  }
  return d;
}

void
SurfaceZPlane::normal(double xyz[3], double uvw[3])
{
  uvw[0] = 0.0;
  uvw[1] = 0.0;
  uvw[2] = 1.0;
}

//==============================================================================

class SurfaceZCylinder : public Surface
{
  double x0, y0, r;
public:
  SurfaceZCylinder(pugi::xml_node surf_node);
  double evaluate(double xyz[3]);
  double distance(double xyz[3], double uvw[3], bool coincident);
  void normal(double xyz[3], double uvw[3]);
};

SurfaceZCylinder::SurfaceZCylinder(pugi::xml_node surf_node)
{
  const char *coeffs = surf_node.attribute("coeffs").value();
  int stat = sscanf(coeffs, "%lf %lf %lf", &x0, &y0, &r);
  if (stat != 3)
  {
    std::cout << "Something went wrong reading surface coeffs!" << std::endl;
  }
}

double
SurfaceZCylinder::evaluate(double xyz[3])
{
  const double x = xyz[0] - x0;
  const double y = xyz[1] - y0;
  double f = x*x + y*y - r*r;
  return f;
}

double
SurfaceZCylinder::distance(double xyz[3], double uvw[3], bool coincident)
{
  double a = 1.0 - uvw[2]*uvw[2];  // u^2 + v^2

  if (a == 0.0) return INFTY;

  double x = xyz[0] - x0;
  double y = xyz[1] - y0;
  double k = x*uvw[0] + y*uvw[1];
  double c = x*x + y*y - r*r;
  double quad = k*k - a*c;

  if (quad < 0.0)
  {
    // No intersection with cylinder.
    return INFTY;
  }
  else if (coincident or fabs(c) < FP_COINCIDENT)
  {
    // Particle is on the cylinder, thus one distance is positive/negative
    // and the other is zero. The sign of k determines if we are facing in or
    // out.
    if (k >= 0.0) return INFTY;
    else return (-k + sqrt(quad)) / a;
  }
  else if (c < 0.0)
  {
    // Particle is inside the cylinder, thus one distance must be negative
    // and one must be positive. The positive distance will be the one with
    // negative sign on sqrt(quad)
    return (-k + sqrt(quad)) / a;
  }
  else
  {
    // Particle is outside the cylinder, thus both distances are either
    // positive or negative. If positive, the smaller distance is the one
    // with positive sign on sqrt(quad)

    double d = (-k - sqrt(quad)) / a;
    if (d < 0.0) return INFTY;
    return d;
  }
}

void
SurfaceZCylinder::normal(double xyz[3], double uvw[3])
{
  uvw[0] = 2.0 * (xyz[0] - x0);
  uvw[1] = 2.0 * (xyz[1] - y0);
  uvw[2] = 0.0;
}

//==============================================================================

extern "C" void
read_surfaces(pugi::xml_node *node)
{
  // Count the number of surfaces.
  int n_surfaces = 0;
  for (pugi::xml_node surf_node = node->child("surface"); surf_node;
       surf_node = surf_node.next_sibling("surface"))
  {
    n_surfaces += 1;
  }

  // Allocate the array of Surface pointers.
  surfaces_c = new Surface* [n_surfaces];

  // Loop over XML surface elements and populate the array.
  {
    pugi::xml_node surf_node;
    int i_surf;
    for (surf_node = node->child("surface"), i_surf = 0; surf_node;
         surf_node = surf_node.next_sibling("surface"), i_surf++)
    {
      if (surf_node.attribute("type"))
      {
        const pugi::char_t *surf_type = surf_node.attribute("type").value();
        if (strcmp(surf_type, "z-cylinder") == 0)
        {
          surfaces_c[i_surf] = new SurfaceZCylinder(surf_node);
        }
        else if (strcmp(surf_type, "z-plane") == 0)
        {
          surfaces_c[i_surf] = new SurfaceZPlane(surf_node);
        }
        else
        {
          std::cout << "Call error or handle uppercase here!" << std::endl;
          std::cout << surf_type << std::endl;
        }
      }
    }
  }
}

//==============================================================================

extern "C" bool
surface_sense(int surf_ind, double xyz[3], double uvw[3])
{
  return surfaces_c[surf_ind]->sense(xyz, uvw);
}

extern "C" double
surface_distance(int surf_ind, double xyz[3], double uvw[3], bool coincident)
{
  return surfaces_c[surf_ind]->distance(xyz, uvw, coincident);
}
