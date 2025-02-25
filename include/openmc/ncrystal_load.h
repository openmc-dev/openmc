#ifndef OPENMC_NCRYSTAL_LOAD_H
#define OPENMC_NCRYSTAL_LOAD_H

#include <functional>
#include <memory>

namespace NCrystalVirtualAPI {

// NOTICE: Do NOT make ANY changes in the NCrystalVirtualAPI::VirtAPI_Type1_v1
// class, it is required to stay exactly constant over time and compatible with
// the same definition used to compile the NCrystal library! But changes to
// white space, comments, and formatting is of course allowed.  This API was
// introduced in NCrystal 4.1.0.

class VirtAPI_Type1_v1 {
public:
  //Note: neutron must be an array of length 4 with values {ekin,ux,uy,uz}
  class ScatterProcess;
  virtual const ScatterProcess* createScatter(const char* cfgstr) const = 0;
  virtual const ScatterProcess* cloneScatter(const ScatterProcess*) const = 0;
  virtual void deallocateScatter(const ScatterProcess*) const = 0;
  virtual double crossSectionUncached(
    const ScatterProcess&, const double* neutron) const = 0;
  virtual void sampleScatterUncached(const ScatterProcess&,
    std::function<double()>& rng, double* neutron) const = 0;
  // Plumbing:
  static constexpr unsigned interface_id = 1001;
  virtual ~VirtAPI_Type1_v1() = default;
  VirtAPI_Type1_v1() = default;
  VirtAPI_Type1_v1(const VirtAPI_Type1_v1&) = delete;
  VirtAPI_Type1_v1& operator=(const VirtAPI_Type1_v1&) = delete;
  VirtAPI_Type1_v1(VirtAPI_Type1_v1&&) = delete;
  VirtAPI_Type1_v1& operator=(VirtAPI_Type1_v1&&) = delete;
};

} // namespace NCrystalVirtualAPI

namespace openmc {

using NCrystalAPI = NCrystalVirtualAPI::VirtAPI_Type1_v1;

std::shared_ptr<const NCrystalAPI> load_ncrystal_api();

class NCrystalScatProc final {
public:
  NCrystalScatProc() {}

  NCrystalScatProc(const char* cfgstr)
    : api_(load_ncrystal_api()), p_(api_->createScatter(cfgstr))
  {}

  // Note: Neutron state array is {ekin,ux,uy,uz}

  double cross_section(const double* neutron_state) const
  {
    return api_->crossSectionUncached(*p_, neutron_state);
  }

  void scatter(std::function<double()>& rng, double* neutron_state) const
  {
    api_->sampleScatterUncached(*p_, rng, neutron_state);
  }

  NCrystalScatProc clone() const
  {
    NCrystalScatProc c;
    if (p_) {
      c.api_ = api_;
      c.p_ = api_->cloneScatter(p_);
    }
    return c;
  }

  // Plumbing (move-only semantics, but supports explicit clone):
  NCrystalScatProc(const NCrystalScatProc&) = delete;
  NCrystalScatProc& operator=(const NCrystalScatProc&) = delete;

  NCrystalScatProc(NCrystalScatProc&& o) : api_(std::move(o.api_)), p_(nullptr)
  {
    std::swap(p_, o.p_);
  }

  NCrystalScatProc& operator=(NCrystalScatProc&& o)
  {
    if (p_) {
      api_->deallocateScatter(p_);
      p_ = nullptr;
    }
    std::swap(api_, o.api_);
    std::swap(p_, o.p_);
    return *this;
  }

  ~NCrystalScatProc()
  {
    if (p_)
      api_->deallocateScatter(p_);
  }

private:
  std::shared_ptr<const NCrystalAPI> api_;
  const NCrystalAPI::ScatterProcess* p_ = nullptr;
};

} // namespace openmc

#endif
