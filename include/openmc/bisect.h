// taken fom boost1.85 - boost::math::tools::bisect
//
// Boost Software License - Version 1.0 - August 17th, 2003

// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:
//
// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#ifndef OPENMC_BISECT_H
#define OPENMC_BISECT_H

#include <iostream>

template <class F, class Tol>
std::pair<double, double> bisect(F f, double min, double max, Tol tol, std::uintmax_t& max_iter)
{
   double fmin = f(min);
   double fmax = f(max);
   if (fmin == 0)
   {
      max_iter = 2;
      return std::make_pair(min, min);
   }
   if (fmax == 0)
   {
      max_iter = 2;
      return std::make_pair(max, max);
   }

   //
   // Error checking:
   //
   if (min >= max)
   {
      throw std::invalid_argument("Arguments min and max in wrong order");
   }
   if (fmin * fmax >= 0)
   {
      throw std::invalid_argument("No change of sign, either there is no root to find, or there are multiple roots in the interval");
   }

   //
   // Three function invocations so far:
   //
   std::uintmax_t count = max_iter;
   if (count < 3) {count = 0;}
   else {count -= 3;}
      

   while (count && (0 == tol(min, max)))
   {
      double mid = (min + max) / 2;
      double fmid = f(mid);
      if ((mid == max) || (mid == min))
         break;
      if (fmid == 0)
      {
         min = max = mid;
         break;
      }
      else if (fmid * fmin < 0)
      {
         max = mid;
      }
      else
      {
         min = mid;
         fmin = fmid;
      }
      --count;
   }

   max_iter -= count;

   return std::make_pair(min, max);
}

#endif