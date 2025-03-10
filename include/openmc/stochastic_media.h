#include <string>
#include <unordered_map>

#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"

#include "openmc/constants.h"
#include "openmc/memory.h" // for unique_ptr
#include "openmc/vector.h"

namespace openmc {
//==============================================================================
// Global variables
//==============================================================================
class Stochastic_Media;

namespace model {
extern std::unordered_map<int32_t, int32_t> stochastic_media_map;
extern vector<unique_ptr<Stochastic_Media>> stochastic_media;
} // namespace model
//==============================================================================
class Stochastic_Media {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors, factory functions

  explicit Stochastic_Media(pugi::xml_node cell_node);
  Stochastic_Media() {};
  virtual ~Stochastic_Media();

  // Accessors
  //! Get name
  //! \return Stochastic Media name
  const std::string& name() const { return name_; }

  //! Set name
  void set_name(const std::string& name) { name_ = name; }

  //! Get ID of stochastic media
  //! \return ID of stochastic media
  int32_t id() const { return id_; }

  //! Assign a unique ID to the stochastic media
  //! \param[in] Unique ID to assign. A value of -1 indicates that an ID
  //!   should be automatically assigned.
  void set_id(int32_t id_);

  //! Get the particle radius
  double radius() const { return radius_; }

  //! Get the packing fraction of the particle material
  double pf() const { return pf_; }

  //! Get the material of the particle
  int32_t particle_mat(int32_t i) const { return particle_mat_[i]; }
  vector<int32_t> particle_mat() const { return particle_mat_; }


  //! Get the material of the matrix
  int32_t matrix_mat() const { return matrix_mat_; }


  //----------------------------------------------------------------------------
  // Data
  int32_t id_ {C_NONE}; //!< Unique ID
  std::string name_;    //!< Name of stochastic_media
  // Parameters: the  particle radius, now this method only support sphere
  double radius_;
  // The packing fraction of the particle material
  double pf_;
  // The material of the particle
  vector<int32_t> particle_mat_;
  // The material of the matrix
  int32_t matrix_mat_;
  //----------------------------------------------------------------------------

  virtual void sample_material(Particle& p) = 0;
  virtual double distance_to_stochamedia(Particle& p) = 0;
  virtual void adjust_indices();

protected:
  // Protected data members
  int64_t index_;
};

class CLS_Media : public Stochastic_Media {
public:
  explicit CLS_Media(pugi::xml_node cell_node);
  CLS_Media() {};

  void sample_material(Particle& p) override;
  double distance_to_stochamedia(Particle& p) override;


};
//==============================================================================
// Non-member functions
//==============================================================================
//! Read stochastic media data XML node
//! \param[in] root node of stochastic_media XML element
void read_stochastic_media(pugi::xml_node root);

void free_memory_stochastic_media();

} // namespace openmc