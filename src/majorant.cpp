#include <fstream>

#include <fmt/core.h>

#include "openmc/constants.h"
#include "openmc/majorant.h"
#include "openmc/material.h"
#include "openmc/nuclide.h"
#include "openmc/search.h"
#include "openmc/simulation.h"
#include "openmc/thermal.h"

namespace openmc {

namespace data {
std::vector<std::unique_ptr<Majorant>> nuclide_majorants;
std::unique_ptr<Majorant> n_majorant;
}

void create_majorant() {
  write_message("Creating majorant cross section...");
  // create a majorant XS for each nuclide
  for (const auto& nuclide : data::nuclides) {
    data::nuclide_majorants.push_back(std::make_unique<Majorant>());
    auto& majorant = data::nuclide_majorants.back();

    for (int t = 0; t < nuclide->kTs_.size(); t++) {
      auto total = nuclide->xs_[t].slice(openmc::tensor::all, 0);
      std::vector<double> xs(total.begin(), total.end());
      auto energies = nuclide->grid_[t].energy;
      majorant->update(energies, xs);

      // include unresolved resonance region data
      // in the majorant if present
      if (nuclide->urr_present_) {
        const auto& urr_data = nuclide->urr_data_[t];
        std::vector<double> energies(urr_data.energy_.begin(), urr_data.energy_.end());

        double max_urr_total {0.0};

        for (auto xs_vals : urr_data.xs_values_) {
          max_urr_total = std::max(max_urr_total, xs_vals.total);
        }

        if(urr_data.interp_ == Interpolation::log_log) {
          std::cout << fmt::format("Nuclide {} uses log-log interpolation", nuclide->name_) << std::endl;
        }

        if (urr_data.multiply_smooth_) {
          std::cout << "Multiplied URR: " << nuclide->name_ << std::endl;
          majorant->grid_.init();
          std::vector<double> xs_vals;
          for (int i = 0; i < energies.size(); i++) {
            xs_vals.push_back(majorant->calculate_xs(energies[i]) * max_urr_total);
          }
          majorant->update(energies, xs_vals);
        } else {
          std::vector<double> xs_vals(energies.size(), max_urr_total);
          majorant->update(energies, xs_vals);
        }
        // majorant->update_urr(energies, xs, urr_data.interp_);
      }
    }
    // initialize the energy grid for this nuclide
    majorant->grid_.init();
    // majorant->write_ascii(nuclide->name_ + "_majorant.txt");
  }

  auto majorant_e_grid = compute_majorant_energy_grid();


  std::vector<double> xs_vals;

  std::vector<Majorant> macro_majorants;

  // compute a majorant for every material
  for (auto& material : model::materials) {
    std::vector<double> material_xs;
    for (auto e_val : majorant_e_grid) {
      double xs_val = 0.0;

      // See if thermal data is present for any of the material nuclides
      bool check_sab = material->thermal_tables_.size() > 0;
      int thermal_table_idx = 0;

      int j = 0;

      for (int i = 0; i < material->nuclide_.size(); i++) {
        int i_nuc = material->nuclide_[i];
        xs_val += material->atom_density_(i) * data::nuclide_majorants[i_nuc]->calculate_xs(e_val);

        int i_sab = C_NONE;
        double sab_frac = 0.0;

        if (check_sab) {
          const auto& sab {material->thermal_tables_[j]};
          if (i == sab.index_nuclide) {
            // Get index in the sab_tables
            i_sab = sab.index_table;
            sab_frac = sab.fraction;

            ++j;
            if (j == material->thermal_tables_.size()) check_sab = false;
          }
        }

        // compute the max sab xs if present
        if (i_sab >= 0) {
          const auto& tdata = data::thermal_scatt[i_sab];
          // compute the maximum xs value over all temperatures
          double thermal_xs = 0.0;
          for (int k = 0; k < tdata->kTs_.size(); k++) {
            // compute the thermal xs at the temperature and energy
            double elastic, inelastic;
            tdata->data_[k].calculate_xs(e_val, &elastic, &inelastic);
            if (elastic + inelastic > thermal_xs) {
              thermal_xs = elastic + inelastic;
            }
          }
          // adjust the current xs_val for the thermal component
          xs_val += sab_frac * thermal_xs;
        }

      }

      material_xs.push_back(xs_val);
    }
    macro_majorants.emplace_back(Majorant());
    macro_majorants.back().update(majorant_e_grid, material_xs);
    macro_majorants.back().write_ascii(fmt::format("mat_{}_majorant.txt", material->id_));
  }

  data::n_majorant = std::make_unique<Majorant>();

  for (auto& macro_majorant : macro_majorants) {
    data::n_majorant->update(macro_majorant.grid_.energy, macro_majorant.xs_);
  }

  data::n_majorant->grid_.init();
  data::n_majorant->write_ascii("macro_majorant.txt");
}

std::vector<double>
compute_majorant_energy_grid() {

  std::vector<double> common_e_grid;
  for (const auto& nuc_maj : data::nuclide_majorants) {
    auto& e_grid = nuc_maj->grid_.energy;
    // append new points to the current group of points
    common_e_grid.insert(common_e_grid.end(), e_grid.begin(), e_grid.end());

    // remove duplicates
    std::unique(common_e_grid.begin(), common_e_grid.end());
  }
  std::sort(common_e_grid.begin(), common_e_grid.end());

  // remove all values below the minimum neutron energy
  int neutron = ParticleType::neutron().transport_index();
  auto min_it = common_e_grid.begin();
  while (*min_it < data::energy_min[neutron]) { min_it++; }
  common_e_grid.erase(common_e_grid.begin(), min_it + 1);
  // insert the minimum neutron energy at the beginning
  common_e_grid.insert(common_e_grid.begin(), data::energy_min[neutron]);

  // remove all values above the maximum neutron energy
  auto max_it = --common_e_grid.end();
  while (*max_it > data::energy_max[neutron]) { max_it--; }
  common_e_grid.erase(max_it - 1, common_e_grid.end());
  // insert the maximum neutron energy at the end
  common_e_grid.insert(common_e_grid.end(), data::energy_max[neutron]);

  return common_e_grid;
}


Majorant::Majorant(const std::vector<double>& energy,
                   const std::vector<double>& xs) : xs_(xs)
{
  grid_.energy = energy;
  grid_.init();
}

double
Majorant::calculate_xs(double energy) const
{
  // Find energy index on energy grid
  int neutron = ParticleType::neutron().transport_index();
  int i_log_union = std::log(energy * data::energy_min_rcp[neutron]) * simulation::log_spacing_rcp;

  int i_grid;
  if (i_log_union < 0) {
    i_grid = 0;
  } else if (i_log_union >= (grid_.grid_index.size() - 2)) {
    i_grid = grid_.energy.size() - 2;
  } else {
    // Determine bounding indices based on which equal log-spaced
    // interval the energy is in
    int i_low  = grid_.grid_index[i_log_union];
    int i_high = grid_.grid_index[i_log_union + 1] + 1;

    // Perform binary search over reduced range
    i_grid = i_low + lower_bound_index(&grid_.energy[i_low], &grid_.energy[i_high], energy);
  }

  // check for rare case where two energy points are the same
  if (grid_.energy[i_grid] == grid_.energy[i_grid + 1]) ++i_grid;

  // calculate interpolation factor
  double f = (energy - grid_.energy[i_grid]) /
              (grid_.energy[i_grid + 1]- grid_.energy[i_grid]);

  double xs = (1.0 - f) * xs_[i_grid] + f * xs_[i_grid + 1];

  return 1.0 * xs;
}

bool Majorant::intersect_2D(std::pair<double, double> p1,
                            std::pair<double, double> p2,
                            std::pair<double, double> p3,
                            std::pair<double, double> p4,
                            std::pair<double, double>& intersection) {

  double denominator = (p4.first - p3.first) * (p1.second - p2.second) -
                        (p1.first - p2.first) * (p4.second - p3.second);

  // if the lines are parallel, return no intersection
  if (fabs(denominator) <= FP_PRECISION) { return false; }

  double numerator = (p3.second - p4.second) * (p1.first - p3.first) +
                     (p4.first - p3.first) * (p1.second - p3.second);

  double t = numerator / denominator;

  if (t < 0.0 || t > 1.0) { return false; }

  // compute intersection location
  double x = p1.first + (p2.first - p1.first) * t;
  double slope = (p2.second - p1.second) / (p2.first - p1.first);
  double y = p1.second + (x - p1.first) * slope;
  intersection = {x, y};
  return true;
}

bool Majorant::is_above(std::pair<double, double> p1,
                        std::pair<double, double> p2,
                        std::pair<double, double> p3) {
  // if the line is vertical, use
  // comparison of x values
  if (fabs(p2.first - p1.first) < FP_PRECISION) {
    return p3.first < p2.first;
  } else {
    double slope = (p2.second - p1.second) / (p2.first - p1.first);
    double val = p1.second + slope * (p3.first - p1.first);
    return val < p3.second;
  }

}

void Majorant::write_ascii(const std::string& filename) const {

  std::ofstream of(filename);

  for (int i = 0; i < xs_.size(); i++) {
    of << grid_.energy[i] << "\t" << xs_[i] << "\n";
  }

  of.close();
}

void Majorant::update(std::vector<double> energy_other,
                      std::vector<double> xs_other) {

  XS xs_a(grid_.energy, xs_);
  XS xs_b(energy_other, xs_other);

  // early exit checks
  if (xs_b.complete()) { return; }

  if (xs_a.complete()) {
    grid_.energy = energy_other;
    xs_ = xs_other;
    return;
  }

  // resulting output of this algorithm
  std::vector<double> e_out, xs_out;
  std::vector<size_t> mask;

  // references we'll use to keep track of which
  // is currently considered the larger cross section
  XS& current_xs = xs_a;
  XS& other_xs = xs_b;

  // if the other cross section starts at a lower energy, its
  // value is considered to be higher
  if (other_xs.get_e() < current_xs.get_e()) {
    std::swap(current_xs, other_xs);
  }

  // if the two cross sections start at the same energy
  // pick the one with the higher xs value
  if (other_xs.get_e() == current_xs.get_e() && other_xs.get_xs() > current_xs.get_xs()) {
    std::swap(current_xs, other_xs);
  }

  // add the first point to the final cross section
  e_out.push_back(current_xs.get_e());
  xs_out.push_back(current_xs.get_xs());
  current_xs++;

  // continue adding points until the other xs min
  // energy is lower than the output minimum energy
  while(current_xs.get_e() < other_xs.get_e() && !current_xs.complete()) {
    e_out.push_back(current_xs.get_e());
    xs_out.push_back(current_xs.get_xs());
    current_xs++;
  }

  // Corner case: first point of the other xs is higher and nearer
  // There is no intersection to find b/c there is no previous point to
  // use for interpolation. Here, we'll insert a point at that energy
  // in the current cross section, insert the point of the other cross
  // section, and swap the two.
  // NOTE: possible problem with computing intersections with vertical slopes here
  if (current_xs.get_e() <= other_xs.get_e() && is_above({e_out.back(), xs_out.back()}, current_xs.get(), other_xs.get())) {
    // insert point on current xs segment
    double slope = (current_xs.get_xs() - xs_out.back()) /
                    (current_xs.get_e() - e_out.back());
    double xs_val = xs_out.back() + slope * (other_xs.get_e() - e_out.back());
    // insert point from other cross section
    e_out.push_back(other_xs.get_e());
    xs_out.push_back(other_xs.get_xs());
    // swap the cross sections and advance
    std::swap(current_xs, other_xs);
    current_xs.advance(e_out.back());
    other_xs.advance(e_out.back());
  }

  // the next point added to the cross section is selected based 4 different cases
  //
  // above: the other cross section point above the line between the last point added and the current xs's next point
  // below: ther other cross section point is below the line betwee the last point added and the current xs's next point
  // nearer: the other cross section point is closer in energy to the last point added
  // farther: the other cross section point is farther in energy from the last point added
  //
  // Case 1: other xs value is above and nearer
  //
  //         Check for an intersection with the line segment between last_pnt and current_next.
  //         If an intersection exists, insert a point in the majorant at that location and
  //         swap the current and other cross sections.
  //
  // Case 2: other xs value is above and farther
  //
  //         Check for an intersection with the line segment between last_pnt and current_next.
  //         Apply the same treatment as Case 1.
  //
  // Case 3: other xs value is below and nearer
  //
  //         Insert a point at the energy of the other xs value on the line segment between
  //         last_pnt and current_next. This point is superfluous and purely for the algorithm
  //         robustness, so the index of this point is added to the mask of indices to be removed
  //         later.
  //
  // Case 4: other xs value is below and farther
  //
  //         Insert the location of current_next to the output cross section and continue

  // start looping over values
  while (true) {

    if (current_xs.complete() || other_xs.complete()) { break; }

    auto current_next = current_xs.get();
    auto other_next = other_xs.get();
    std::pair<double, double> last_pnt = {e_out.back(), xs_out.back()};

    bool above = is_above(last_pnt, current_next, other_next);
    bool nearer = other_next.first < current_next.first;

    // for clarity
    bool below = !above;
    bool farther = !nearer;

    // Cases 1 & 2
    if (above) {
      // setup intersection check
      auto p1 = last_pnt;
      auto p2 = current_next;
      auto p3 = other_xs.prev();
      auto p4 = other_next;

      std::pair<double, double> intersection;
      auto intersection_found = intersect_2D(p1, p2, p3, p4, intersection);

      if (intersection_found) {
        e_out.push_back(intersection.first);
        xs_out.push_back(intersection.second);
        std::swap(current_xs, other_xs);
      } else {
        e_out.push_back(current_next.first);
        xs_out.push_back(current_next.second);
      }
    // Case 3
    } else if (below and nearer) {
      // compute the y value
      double slope = (current_next.second - last_pnt.second) /
                      (current_next.first - last_pnt.first);
      double xs_val = last_pnt.second + slope * (other_next.first - last_pnt.first);
      e_out.push_back(other_next.first);
      xs_out.push_back(xs_val);
      // add index of this point to the mask
      mask.push_back(e_out.size() - 1);
    // Case 4
    } else if (below and farther) {
      e_out.push_back(current_next.first);
      xs_out.push_back(current_next.second);
    }

    current_xs.advance(e_out.back());
    other_xs.advance(e_out.back());
  }

  // remove any superfluous points
  std::sort(mask.begin(), mask.end(), std::greater<double>());
  for (auto idx : mask) {
    e_out.erase(e_out.begin() + idx);
    xs_out.erase(xs_out.begin() + idx);
  }

  // add any additional data to the output xs
  while (!xs_a.complete()) {
    e_out.push_back(xs_a.get_e());
    xs_out.push_back(xs_a.get_xs());
    xs_a++;
  }

  while (!xs_b.complete()) {
    e_out.push_back(xs_b.get_e());
    xs_out.push_back(xs_b.get_xs());
    xs_b++;
  }

  // increase all xs values in the majorant by 5 percent
  std::transform(xs_out.begin(), xs_out.end(), xs_out.begin(), [](double xs_val) { return xs_val * 1.05; });

  // update the values of the majorant
  grid_.energy = std::move(e_out);
  xs_ = std::move(xs_out);

}

// Majorant XS definitions

std::pair<double, double>
Majorant::XS::get() const { return {energies_.at(idx_), total_xs_.at(idx_)}; }

double Majorant::XS::get_e() const { return energies_.at(idx_); }

double Majorant::XS::get_xs() const { return total_xs_.at(idx_); }

std::pair<double, double>
Majorant::XS::prev() const { return {energies_.at(idx_ - 1), total_xs_.at(idx_ - 1)}; }

double Majorant::XS::prev_e() const { return energies_.at(idx_ - 1); }

double Majorant::XS::prev_xs() const { return total_xs_.at(idx_ - 1); }

void Majorant::XS::advance(double energy) {
  double e = energies_[idx_];
  while (e <= energy && !this->complete()) { e = energies_[++idx_]; }
}

bool Majorant::XS::complete() const { return idx_ >= energies_.size(); }

}
