.. _dnp_drift:

Delayed Neutron Precursor (DNP) Drift
=====================================

.. note::
   The DNP drift feature in OpenMC is currently limited to regular meshes
   (:class:`openmc.RegularMesh`).

.. warning::
   Because DNP removal is an analog process, it is not taken into account with
   collision or track-length processes. One way to correct the averaged
   :math:`k_{\text{eff}}` is to normalize it by the ratio of analog and track-length
   nu-fission tallies.

In solid-fueled reactors, delayed neutron precursors (DNPs) are produced and decay at the
same spatial location where the fission event occurred. In liquid-fueled reactors (e.g.,
molten salt reactors (MSRs)), the fissile material is dissolved in a carrier fluid that
circulates through the primary loop. As a result, DNPs are carried away from their point
of origin by the flowing fuel and may emit their delayed neutrons outside of the active
core region. The spatial redistribution of precursors alters the delayed neutron source
distribution, which can significantly affect reactor criticality and kinetics. This
phenomenon is commonly referred to as **DNP drift**. 

OpenMC provides capabilities to model DNP drift by explicitly transporting precursors along
fluid velocity streamlines. This feature is designed to support high-fidelity multiphysics
coupling with computational fluid dynamics (CFD) codes such as 
`nekRS <https://github.com/Nek5000/nekrs>`_ through the
`Cardinal <https://cardinal.cels.anl.gov/>`_ framework.

Streamline integration
----------------------

In OpenMC, individual delayed neutron precursors are advected along fluid streamlines
until they reach their sampled decay time. The decay time for a precursor belonging to
delayed group :math:`g` is sampled as:

.. math::

   t_{\text{decay}} = - \frac{\log(\zeta)}{\lambda_g},

where :math:`\zeta` is a randomly generated number and :math:`\lambda_g` is the decay
constant corresponding to the sampled delayed neutron group number :math:`g`.

Each precursor trajectory is obtained by integrating:

.. math::

   \frac{dy}{dt} = v(y,t),

where :math:`y(t)` is the precursor position and :math:`v(y,t)` is the local fluid velocity
interpolated from the underlying mesh. In the current implementation, the velocity field
is assumed to be static within an OpenMC calculation, so the time dependence is dropped
and :math:`v(y,t)` reduces to :math:`v(y)`. The integration is performed using the classical
fourth-order Runge-Kutta method (RK4). Given a precursor at position :math:`x_n` and
a time step :math:`\Delta t`, the method evaluates four intermediate velocity estimates:

.. math::

   k_1 = v(y_n), \\
   k_2 = v\left(y_n + \frac{\Delta t}{2}\,k_1\right), \\
   k_3 = v\left(y_n + \frac{\Delta t}{2}\,k_2\right), \\
   k_4 = v\left(y_n + \Delta t\,k_3\right).

The position is then advanced as:

.. math::

   y_{n+1} = y_n + \frac{\Delta t}{6}\left(k_1 + 2\,k_2 + 2\,k_3 + k_4\right).

Boundary conditions and recycling
---------------------------------

In practice, the velocity mesh used for DNP transport typically covers only a portion
of the full primary loop (e.g., the active core region). Some surfaces of the mesh
therefore represent the inlet and outlet of the modeled open system.
These boundary faces are organized into physical groups, following conventions used
by mesh generators and CFD codes. The user must specify which physical groups correspond
to **inlet**, **outlet**, or **wall** boundaries. Additionally, because the DNP drift feature
is currently limited to OpenMC regular meshes (:class:`openmc.RegularMesh`), the assignment
of individual faces to physical groups must also be provided explicitly by the user.

During streamline integration, a precursor may encounter one of these boundaries.
Each type is handled as follows:

- **Outlet**: The precursor exits the modeled domain and is discarded.

- **Wall**: The precursor is stopped at the position where the wall crossing would have
  occurred and decays at that location.

- **Inlet**: The precursor is repositioned at a location sampled uniformly on the set
  of inlet faces and continues its trajectory.

In a real MSR, however, the fuel salt recirculates: it exits the core, passes through
heat exchangers and piping, and re-enters the core at the inlet. Precursors leaving
the core might re-enter if they survive long enough. The **recycling** option
models this recirculation without requiring the velocity mesh to cover the entire
primary loop. It assumes a fixed external travel time :math:`\tau_{\text{ext}}`
representing the transit duration through the portion of the loop that is not
explicitly represented with the velocity mesh.

When recycling is enabled, the algorithmic behavior associated with the outlet boundary
condition is modified: instead of being immediately discarded, a precursor that reaches
an outlet has its remaining lifetime :math:`t_{\text{rem}}` compared to
:math:`\tau_{\text{ext}}`:

- If :math:`t_{\text{rem}} \gt \tau_{\text{ext}}`, the precursor survives the
  external loop. It is re-injected at a randomly sampled inlet position with an updated
  remaining lifetime of :math:`t_{\text{rem}} - \tau_{\text{ext}}`, and streamline
  integration resumes from the new inlet location.

- If :math:`t_{\text{rem}} \le \tau_{\text{ext}}`, the precursor decays during
  the external transit and is discarded.

This cycle of exit, external transit, and re-injection may repeat multiple times for
long-lived precursors, allowing them to complete several passes through the core
before eventually decaying.

Geometric reconciliation
------------------------

After transport of each DNP, a **reconciliation step** verifies that the
precursor site lies within the OpenMC geometry. Because DNP transport is
performed on a dedicated mesh while particle transport uses the native OpenMC
geometry (CSG or DAGMC), a precursor may be transported to a position that the
geometry routines classify as outside the model.

Note that this inside/outside classification depends on the **neutron emission direction**
sampled during the original fission event, not the fluid velocity vector. A site physically
located at a boundary surface may be seen as inside or outside the OpenMC geometry depending
on whether this direction points inward or outward relative to the surface normal.

The reconciliation procedure distinguishes two cases:

1. **Site classified as inside the geometry** - OpenMC's cell-finding routines successfully
   locate the precursor inside a cell (including sites on a boundary surface whose
   emission direction points inward). The site is accepted as a valid delayed neutron source
   location for the next generation.

2. **Site not found inside any cell** - The precursor is nudged slightly backward
   along its emission direction and the cell search is repeated:

   - **Site on a boundary surface (heading outward)**: If the nudged position is
     found inside a cell and the nearest boundary is within the nudge tolerance, 
     the precursor is determined to be on a surface. The boundary condition associated
     with that surface is then applied:

      - *Reflective / White*: the site position and direction are transformed according
        to the boundary condition. If the boundary carries an albedo factor, the site
        weight is adjusted accordingly.
      - *Periodic*: the site is mapped to the corresponding periodic partner surface.
      - *Vacuum*: the site is discarded.

   - **Site genuinely outside**: If the nudged position is still outside the geometry,
     or the nearest boundary is closer than the nudge tolerance, the site is
     discarded. A warning is emitted so that the user can diagnose potential
     issues with mesh alignment.

.. tip::
   A high number of lost sites during reconciliation might indicate a mismatch between
   the velocity mesh and the OpenMC geometry, or an integrator time step that is too large.
   Reducing ``integrator_dt`` or improving the geometric alignment between the mesh and
   the model can mitigate this.

Usage
-----

The following example shows how to set up DNP drift on a simple 2x2x2 regular mesh.

**1. Define the mesh and the velocity field**

.. code-block:: python

   import numpy as np
   import openmc

   # Create a 2x2x2 regular mesh
   mesh = openmc.RegularMesh()
   mesh.lower_left = (0.0, 0.0, 0.0)
   mesh.upper_right = (10.0, 10.0, 10.0)
   mesh.dimension = (2, 2, 2)

   # Define a uniform velocity field on the mesh nodes
   # For a 2x2x2 mesh there are 3x3x3 = 27 unique nodes, each with 3 components.
   velocity_field = openmc.VelocityField(
      mesh=mesh,
      values=np.array([1.0]*3*27).reshape(-1,3),
      mapping="nodal"
   )

**2. Create the DNPDrift object**

.. code-block:: python

   dnp_drift = openmc.DNPDrift(
      velocity_field=velocity_field,
      boundary_map={
         "inlet": [1],
         "outlet": [2],
         "wall": [3]
      },
      physical_group_map={
         "face_ids": [
               0, 24, 12, 36, 7, 31, 19, 43, 15, 39, 21, 45, 2,
               26, 8,32, 4, 16, 10, 22, 29, 41, 35, 47
         ],
         "physical_groups": [
               1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
               3, 3, 3, 3, 3, 3, 3, 3
         ]},
      integrator="RK4",
      integrator_dt=0.1,
      recycling=True,
      external_travel_time=2.0
   )

The parameters are explained below:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Parameter
     - Description
   * - ``velocity_field``
     - The velocity field used to transport precursor sites.
   * - ``boundary_map``
     - Dictionary mapping boundary type names (``inlet``, ``outlet``,
       ``wall``) to lists of physical group numbers.
   * - ``physical_group_map``
     - Dictionary with two keys: ``face_ids`` (integer IDs of each boundary
       face on the mesh) and ``physical_groups`` (the corresponding physical
       group number for each face).
   * - ``integrator``
     - Time integration scheme. Currently, only ``"RK4"`` is supported.
   * - ``integrator_dt``
     - Time step size (in seconds) used by the integrator.
   * - ``recycling``
     - If ``True``, precursors that exit through an outlet can be re-injected at
       an inlet if they do not decay during the external travel. Otherwise, precursors
       leaving by an outlet are considered lost.
   * - ``external_travel_time``
     - Average time (in seconds) needed for a precursor to re-enter the modeled system
       after having left the system by crossing an outlet. Only used when ``recycling=True``. 

**3. Register the DNP drift model in the settings**

.. code-block:: python

   settings = openmc.Settings()
   settings.dnp_drift = dnp_drift
