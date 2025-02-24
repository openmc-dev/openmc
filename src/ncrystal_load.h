#ifndef ncrystal_load
#define ncrystal_load

#include <functional>
#include <memory>

namespace NCrystalDynamicAPI {

  class DynAPI_Type1_v1 {
  public:

    //NOTICE: Do NOT make ANY changes in the DynAPI_Type1_v1 class, it is
    //required to stay exactly constant over time and compatible with the same
    //definition used to compile the NCrystal library! But changes to white
    //space, comments, and formatting is of course allowed.
    //This API was introduced in NCrystal 4.1.0.

    static constexpr unsigned interface_id = 1001;//1000*typenumber+version

    class ScatterProcess;
    virtual const ScatterProcess * createScatter( const char * cfgstr ) const = 0;
    virtual void deallocateScatter( const ScatterProcess * ) const = 0;

    //NB: Cross section units returned are barn/atom:
    virtual double crossSectionUncached( const ScatterProcess&,
                                         double neutron_ekin_eV,
                                         double neutron_dir_ux,
                                         double neutron_dir_uy,
                                         double neutron_dir_uz ) const = 0;
    virtual void sampleScatterUncached( const ScatterProcess&,
                                        std::function<double()>& rng,
                                        double& neutron_ekin_eV,
                                        double& neutron_dir_ux,
                                        double& neutron_dir_uy,
                                        double& neutron_dir_uz ) const = 0;

    virtual ~DynAPI_Type1_v1() = default;
    DynAPI_Type1_v1() = default;
    DynAPI_Type1_v1( const DynAPI_Type1_v1& ) = delete;
    DynAPI_Type1_v1& operator=( const DynAPI_Type1_v1& ) = delete;
    DynAPI_Type1_v1( DynAPI_Type1_v1&& ) = delete;
    DynAPI_Type1_v1& operator=( DynAPI_Type1_v1&& ) = delete;
  };

}

namespace openmc {

  using NCrystalAPI = NCrystalDynamicAPI::DynAPI_Type1_v1;

  std::shared_ptr<const NCrystalAPI> load_ncrystal_api();

  class NCrystalScatProc final {
  public:
  NCrystalScatProc( const char * cfgstr )
    : m_api( load_ncrystal_api() ),
      m_p( m_api->createScatter( cfgstr ) )
      {
      }
    ~NCrystalScatProc()
      {
        if ( m_p )
          m_api->deallocateScatter( m_p );
      }

    struct NeutronState { double ekin, ux, uy, uz; };

    double cross_section( const NeutronState& n ) const
    {
      return m_api->crossSectionUncached( *m_p, n.ekin, n.ux, n.uy, n.uz );
    }
    void scatter( std::function<double()>& rng, NeutronState& n ) const
    {
      m_api->sampleScatterUncached( *m_p, rng, n.ekin, n.ux, n.uy, n.uz );
    }

    NCrystalScatProc( const NCrystalScatProc& ) = delete;
    NCrystalScatProc& operator=( const NCrystalScatProc& ) = delete;
  NCrystalScatProc( NCrystalScatProc&& o  )
    : m_api(std::move(o.m_api)), m_p(nullptr)
      {
        std::swap( m_p, o.m_p );
      }
    NCrystalScatProc& operator=( NCrystalScatProc&& o )
      {
        if ( m_p ) {
          m_api->deallocateScatter( m_p );
          m_p = nullptr;
        }
        std::swap( m_api, o.m_api );
        std::swap( m_p, o.m_p );
        return *this;
      }

  private:
    std::shared_ptr<const NCrystalAPI> m_api;
    const NCrystalAPI::ScatterProcess * m_p;
  };

}

#endif
