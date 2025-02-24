#include "ncrystal_load.hh"
#include <string>
#include <mutex>
#include <stdexcept>
#include <stdio.h>

#if !defined(NCLOAD_WINDOWS) && ( defined (_WIN32) || defined (WIN32) )
#  define NCLOAD_WINDOWS
#endif
#ifdef NCLOAD_WINDOWS
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#else
#  include <dlfcn.h>
#endif

namespace openmc {
  namespace {

    struct NCrystalConfig {
      std::string shlibpath;
      unsigned long intversion = 0;
      std::string symbol_namespace;
    };

    NCrystalConfig query_ncrystal_config()
    {
      char buffer[4096];
#ifdef NCLOAD_WINDOWS
      FILE* pipe = _popen("ncrystal-config --show "
                          "intversion shlibpath namespace", "r");
#else
      FILE* pipe = popen("ncrystal-config --show "
                         "intversion shlibpath namespace 2>/dev/null", "r");
#endif
      if (!pipe)
        return {};//failure
      auto readLine = [pipe]( std::string& tgt ) -> bool
      {
        //Read line and discard trailing whitespace (including newline chars).
        char buffer[4096];
        if (fgets(buffer, sizeof(buffer), pipe) == NULL)
          return false;
        tgt = buffer;
        while( !tgt.empty() && std::isspace( tgt.back() ) )
          tgt.pop_back();
        return true;
      };
      auto parseIntVersion = []( const std::string& s )
      {
        char * str_end = nullptr;
        unsigned long v = std::strtoul( s.c_str(), &str_end, 10 );
        return ( v >= 2002000
                 && v < 999999999
                 && str_end == s.c_str() + s.size() ) ? v : 0;
      };

      NCrystalConfig res;
      bool all_ok(true);
      if ( !readLine(res.shlibpath)
           || !(res.intversion = parseIntVersion( res.shlibpath ))
           || !readLine( res.shlibpath )
           || res.shlibpath.empty()
           || !readLine(res.symbol_namespace ) ) {
        res.intversion = 0;//failure
      }

#ifdef NCLOAD_WINDOWS
      auto returnCode = _pclose(pipe);
#else
      auto returnCode = pclose(pipe);
#endif
      if ( returnCode == 0 && res.intversion >= 2002000  )
        return res;
      return {};//failure
    }

    struct NCrystalAPIDB {
      std::mutex mtx;
      std::shared_ptr<const NCrystalAPI> api;
      typedef void* (*FctSignature)(int);
      FctSignature ncrystal_access_dynapi_fct = nullptr;
    };

    void* load_dynapi_raw( unsigned interface_id, NCrystalAPIDB& db )
    {
      if ( !db.ncrystal_access_dynapi_fct ) {
        auto cfg = query_ncrystal_config();
        if ( !cfg.intversion >= 4000003 )
          throw std::runtime_error("Could not locate a functioning and"
                                   " recent enough NCrystal installation.");
#ifdef NCLOAD_WINDOWS
        auto handle = LoadLibrary(cfg.shlibpath.c_str());
#else
        dlerror();//clear previous errors
        void * handle = dlopen( cfg.shlibpath.c_str(), RTLD_LOCAL | RTLD_LAZY );
#endif
        if (!handle)
          throw std::runtime_error("Loading of the NCrystal library failed");

        std::string symbol("ncrystal");
        symbol += cfg.symbol_namespace;
        symbol += "_access_dynamic_api";

#ifdef NCLOAD_WINDOWS
        FARPROC fproc;
        void * addr = (void *)(intptr_t)GetProcAddress(handle,
                                                       symbol.c_str());
        if (!addr)
          throw std::runtime_error("GetProcAddress("
                                   "ncrystal_access_dynamic_api) failed");
#else
        dlerror();//clear previous errors
        void * addr = dlsym( handle, symbol.c_str());
        if ( !addr )
          throw std::runtime_error("dlsym(ncrystal_access_dynamic_api) failed");
#endif
        db.ncrystal_access_dynapi_fct
          = reinterpret_cast<NCrystalAPIDB::FctSignature>(addr);
      }

      return (*db.ncrystal_access_dynapi_fct)( interface_id );
    }

    NCrystalAPIDB& get_ncrystal_api_db()
    {
      static NCrystalAPIDB db;
      return db;
    }
  }

  std::shared_ptr<const NCrystalAPI> load_ncrystal_api()
  {
    auto& db = get_ncrystal_api_db();
    std::lock_guard<std::mutex> lock(db.mtx);
    if ( ! db.api ) {
      void * raw_api = load_dynapi_raw(NCrystalAPI::interface_id,db);
      if (!raw_api)
        throw std::runtime_error("Failed to access required NCrystal interface.");
      db.api = *reinterpret_cast<std::shared_ptr<const NCrystalAPI>*>(raw_api);
    }
    return db.api;
  }
}
