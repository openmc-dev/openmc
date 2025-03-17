#ifndef OPENMC_TALLIES_FILTER_WEIGHTS_H
#define OPENMC_TALLIES_FILTER_WEIGHTS_H

#include "openmc/tallies/filter.h"
#include "openmc/vector.h"
#include "openmc/span.h"

namespace openmc {

//==============================================================================
//! Bins the weights of the particles.
//==============================================================================

class WeightFilter : public Filter{
public: 
    //----------------------------------------------------------------------------
    // Constructors, destructors
    
    ~WeightFilter() = default;
    
    //----------------------------------------------------------------------------
    // Methods
    
    std::string type_str() const override { return "weight"; }
    FilterType type() const override { return FilterType::WEIGHT; }
    
    void from_xml(pugi::xml_node node) override;

    void get_all_bins(const Particle& p, TallyEstimator estimator,
        FilterMatch& match) const override;
    
    void to_statepoint(hid_t filter_group) const override;
    
    std::string text_label(int bin) const override;
    
    //----------------------------------------------------------------------------
    // Accessors

    const vector<double>& bins() const { return bins_; }
    void set_bins(span<const double> bins);

protected:
    //----------------------------------------------------------------------------
    // Data members
    vector<double> bins_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_WEIGHTS_H