#ifndef XTL_JSON_HPP
#define XTL_JSON_HPP

#include "nlohmann/json.hpp"
#include "xoptional.hpp"

namespace xtl
{
    /***********************************************************
     * to_json and from_json specialization for xtl::xoptional *
     ***********************************************************/

    template <class D>
    void to_json(nlohmann::json& j, const xoptional<D>& o)
    {
        if (!o.has_value())
        {
            j = nullptr;
        }
        else
        {
            j = o.value();
        }
    }

    template <class D>
    void from_json(const nlohmann::json& j, xoptional<D>& o)
    {
        if (j.is_null())
        {
            o = missing<D>();
        }
        else
        {
            o = j.get<D>();
        }
    }
}

#endif
