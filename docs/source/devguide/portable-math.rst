.. _devguide_portable_math:

================================
Reproducible Transcendental Math
================================

OpenMC's regression tests compare simulation output against reference results
bit for bit. That only works when every platform computes the same bits. The
basic IEEE 754 operations (add, subtract, multiply, divide, square root and
fused multiply-add) are correctly rounded everywhere, but transcendental
functions such as ``exp``, ``log``, ``sin``, ``cos`` and ``pow`` are provided by
the system math library, and glibc, Apple's libSystem and the Microsoft C
runtime each implement them differently. Even two glibc versions differ. A
difference in the last bit of a single ``log`` call changes a sampled collision
distance, which changes whether a particle crosses a surface, and from there
the two histories diverge completely. This is why the regression tests fail on
macOS and Windows, and on Linux systems with a newer glibc than the one the
references were generated with (see issue
`#3820 <https://github.com/openmc-dev/openmc/issues/3820>`_).

OpenMC addresses this in two parts.

1. The C++ transport code calls transcendental functions through the ``openmc``
   namespace, for example ``openmc::exp`` or ``openmc::pow``, by including
   ``openmc/math.h``. By default these are the standard library functions, so
   the default build is unchanged.

2. Configuring with ``-DOPENMC_USE_COREMATH=on`` maps the same names to the
   `CORE-MATH <https://core-math.gitlabpages.inria.fr/>`_ project's
   implementations. CORE-MATH provides *correctly rounded* versions of the
   standard math functions: for every input, the result is the double closest
   to the exact mathematical value. Correct rounding fixes the result uniquely,
   so every platform and compiler that builds with it produces the same bits by
   construction for these functions. System libraries are often, but not
   always, correctly rounded; the table below shows where glibc differs on the
   sampled inputs.
   The CI reproducibility workflow builds this way on x86-64 Linux, arm64 Linux
   and arm64 macOS and checks that the regression tests and a digest of every
   particle state of a simulation agree across the three.

Rules for Transport Code
------------------------

- Call transcendental functions through the ``openmc`` namespace and include
  ``openmc/math.h``. Never call ``std::exp``, ``std::log``, ``std::pow``,
  ``std::sin``, ``std::cos``, ``std::tan``, ``std::asin``, ``std::acos``,
  ``std::atan``, ``std::atan2``, ``std::sinh``, ``std::cosh``, ``std::erf``,
  ``std::erfc``, ``std::lgamma``, ``std::tgamma``, ``std::expm1``,
  ``std::log1p``, ``std::log2``, ``std::cbrt`` or ``std::hypot`` from code
  under ``src`` or ``include/openmc``. The unit test
  ``tests/unit_tests/test_portable_math_usage.py`` fails when such a call is
  found, so a new call site cannot silently reintroduce a dependence on the
  system library.
- ``std::sqrt``, ``std::fabs``, ``std::floor``, ``std::ceil``, ``std::fmod``,
  ``std::fma`` and the other basic operations are correctly rounded by the
  IEEE 754 standard and may be used directly.
- Three places still reach the system library and are outside the guarantee:
  the Faddeeva function under ``src/external`` (Doppler broadening with
  windowed multipole data), the quartic solver under ``src/external`` (torus
  surfaces), and the Bessel function ``std::cyl_bessel_j`` used by the
  parametric tokamak source. Code under ``src/external`` is bundled from other
  projects and is exempt from the rule above.

The CORE-MATH Backend
---------------------

CORE-MATH has no releases and no build system; each function is one
self-contained C file. ``cmake/CoreMath.cmake`` fetches the files OpenMC
needs (21 functions and their headers) from the upstream repository at a
pinned commit, verifies each against a recorded SHA-256 hash, and compiles
them into a static library linked into ``libopenmc``. No CORE-MATH code is
kept in the OpenMC repository; updating the backend means changing the commit
and the hashes in that file and regenerating the known-answer vectors
described below.

Build options:

``OPENMC_USE_COREMATH``
  Enables the backend. Requires GCC or Clang (CORE-MATH uses 128-bit integer
  arithmetic, which MSVC does not provide).

``OPENMC_COREMATH_SOURCE_DIR``
  A local checkout of the CORE-MATH repository to compile from instead of
  downloading, for builds without network access. The files must match the
  pinned commit; their hashes are checked the same way as for downloads.

``OPENMC_COREMATH_FMA``
  On x86-64, whether to compile CORE-MATH with the FMA instruction (``AUTO``,
  ``ON`` or ``OFF``). CORE-MATH calls ``fma`` heavily, and on x86-64 the
  baseline instruction set has no FMA, so without the instruction each call
  goes to a software routine and the functions are several times slower. The
  results are identical either way, because ``fma`` is correctly rounded in
  both cases. ``AUTO`` (the default) uses the instruction only when the
  compiler flags already target it, for example with ``-march=native``, so
  the binary never requires more than the rest of the build does; ``ON``
  forces it and the binary then needs an FMA-capable CPU.

The sources are compiled with FMA contraction disabled, because CORE-MATH's
error-free transformations assume it is off; its explicit ``fma`` calls are
unaffected.

Verification
------------

Correct rounding is a property that can be checked, not just claimed. Two
tools in ``tools/dev`` do so:

- ``portable_math_sample.cpp`` evaluates the ``openmc::`` functions and the
  system library on inputs sampled uniformly over the IEEE 754 bit patterns
  of each function's domain, together with a fixed set of special values, and
  records every result as raw bits.
- ``portable_math_report.py`` reads those samples, computes the correctly
  rounded reference for each input with `gmpy2 <https://pypi.org/project/gmpy2/>`_
  (MPFR), and reports the error in units in the last place (ulp) for both
  implementations.

To reproduce the report from a build directory configured with
``OPENMC_USE_COREMATH=on`` (the ``coremath`` static library is in
``<build>/lib/libcoremath.a``)::

    g++ -O2 -std=c++17 -DOPENMC_COREMATH -Iinclude \
        tools/dev/portable_math_sample.cpp <build>/lib/libcoremath.a \
        -o portable_math_sample
    ./portable_math_sample --points 1000000 | \
        python tools/dev/portable_math_report.py - --markdown report.md

The table below is the result for the pinned CORE-MATH commit, one million
bit-uniform samples per function plus the special values. The system library
column is glibc on x86-64, for comparison.

.. include:: portable-math-accuracy.txt

Testing
-------

``tests/cpp_unit_tests/test_portable_math.cpp`` checks the backend against
known-answer vectors in ``portable_math_kat.h``, which record the correctly
rounded bits for the special values and a sample of ordinary inputs of all 21
functions. The test fails on any platform or compiler that produces a
different bit, so a reproducibility problem is caught by ``ctest`` before it
can appear as a regression test difference. The test is skipped in builds
without the backend. The vectors are regenerated from a sample file when the
pinned commit changes::

    python tools/dev/portable_math_report.py samples.txt \
        --kat tests/cpp_unit_tests/portable_math_kat.h
