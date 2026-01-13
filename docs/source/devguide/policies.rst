.. _devguide_policies:

========
Policies
========

---------------------
Python Version Policy
---------------------

OpenMC follows the Scientific Python Ecosystem Coordination guidelines `SPEC 0
<https://scientific-python.org/specs/spec-0000/>`_ on minimum supported
versions, which recommends that support for Python versions be dropped 3 years
after their initial release.

-------------------
C++ Standard Policy
-------------------

C++ code in OpenMC must conform to the most recent C++ standard that is fully
supported in the `version of the gcc compiler
<https://gcc.gnu.org/projects/cxx-status.html>`_ that is distributed with the
oldest version of Ubuntu that is still within its `standard support period
<https://ubuntu.com/about/release-cycle>`_. Ubuntu 22.04 LTS will be supported
through April 2027 and is distributed with gcc 11.4.0, which fully supports the
C++17 standard.

--------------------
CMake Version Policy
--------------------

Similar to the C++ standard policy, the minimum supported version of CMake
corresponds to whatever version is distributed with the oldest version of Ubuntu
still within its standard support period. Ubuntu 22.04 LTS is distributed with
CMake 3.22.
