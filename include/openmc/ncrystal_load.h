#ifndef ncrystal_load
#define ncrystal_load

#include <functional>
#include <memory>

namespace NCrystalDynamicAPI {

// NOTICE: Do NOT make ANY changes in the NCrystalDynamicAPI::DynAPI_Type1_v1
// class, it is required to stay exactly constant over time and compatible
// with the same definition used to compile the NCrystal library! But changes
// to white space, comments, and formatting is of course allowed.  This API
// was introduced in NCrystal 4.1.0.

class DynAPI_Type1_v1 {
public:
  static constexpr unsigned interface_id = 1001; // 1000*typenumber+version

  class ScatterProcess;
  virtual const ScatterProcess* createScatter(const char* cfgstr) const = 0;
  virtual void deallocateScatter(const ScatterProcess*) const = 0;

  // NB: Cross section units returned are barn/atom:
  virtual double crossSectionUncached(const ScatterProcess&,
    double neutron_ekin_eV, double neutron_dir_ux, double neutron_dir_uy,
    double neutron_dir_uz) const = 0;
  virtual void sampleScatterUncached(const ScatterProcess&,
    std::function<double()>& rng, double& neutron_ekin_eV,
    double& neutron_dir_ux, double& neutron_dir_uy,
    double& neutron_dir_uz) const = 0;

  virtual ~DynAPI_Type1_v1() = default;
  DynAPI_Type1_v1() = default;
  DynAPI_Type1_v1(const DynAPI_Type1_v1&) = delete;
  DynAPI_Type1_v1& operator=(const DynAPI_Type1_v1&) = delete;
  DynAPI_Type1_v1(DynAPI_Type1_v1&&) = delete;
  DynAPI_Type1_v1& operator=(DynAPI_Type1_v1&&) = delete;
};

} // namespace NCrystalDynamicAPI

namespace openmc {

using NCrystalAPI = NCrystalDynamicAPI::DynAPI_Type1_v1;

std::shared_ptr<const NCrystalAPI> load_ncrystal_api();

class NCrystalScatProc final {
public:
  NCrystalScatProc(const char* cfgstr)
    : api_(load_ncrystal_api()), p_(api_->createScatter(cfgstr))
  {
  }

  struct NeutronState { double ekin, ux, uy, uz; };

  double cross_section(const NeutronState& n) const
  {
    return api_->crossSectionUncached(*p_, n.ekin, n.ux, n.uy, n.uz);
  }

  void scatter(std::function<double()>& rng, NeutronState& n) const
  {
    api_->sampleScatterUncached(*p_, rng, n.ekin, n.ux, n.uy, n.uz);
  }

  //Plumbing (move-only semantics):
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
  const NCrystalAPI::ScatterProcess* p_;
};

} // namespace openmc

#endif
