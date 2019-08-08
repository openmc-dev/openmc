#ifndef MEASURE_HPP

#include "moab/CN.hpp"

double edge_length( const double* start_vtx_coords,
                    const double* end_vtx_coords );

double measure( moab::EntityType type,
                int num_vertices,
                const double* vertex_coordinatee );

#endif

                
