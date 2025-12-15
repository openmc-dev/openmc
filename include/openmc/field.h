#ifndef OPENMC_FIELD_H
#define OPENMC_FIELD_H

#include "openmc/vector.h"
#include "openmc/mesh.h"

namespace openmc {

class TemperatureField{

public:
    Mesh* mesh_ptr;
    vector<double> values;

    TemperatureField(){};

    TemperatureField(Mesh* mesh_ptr, vector<double> values);

    double distance_to_next_cell(Position r, Direction d);

    double get_temperature(Position r);

    double get_sqrtkT(Position r);

    // GeometryState, particle attribute: temperature mesh bin

};

} //  namespace openmc
#endif // OPENMC_FIELD_H
