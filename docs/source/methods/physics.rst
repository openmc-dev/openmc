.. _methods_physics:

=======
Physics
=======

------------------------------------------
Effect of Thermal Motion on Cross-Sections
------------------------------------------

When a neutron scatters off of a nucleus, many times it is assumed that the
target nucleus is at rest. However, if the material is at a temperature greater
than 0 K, it will have motion associated with the thermal vibration. Thus, the
velocity of the neutron relative to the target nucleus is in general not the
same as the velocity of the neutron entering the collision.

The affect of the thermal motion on the interaction probability can be written
as

.. math::
    :label: freegas1

    v_n \sigma (v_n, T) = \int_0^\infty d\mathbf{v}_T \sigma(v_r, 0)
    \mathbf{v}_r p(\mathbf{v}_T)
    
where :math:`v_n` is the magnitude of the velocity of the neutron,
:math:`\mathbf{v}_T` is the velocity of the target nucleus, :math:`\mathbf{v}_r`
is the relative velocity, and :math:`T` is the temperature of the target
material. In a Monte Carlo code, one must account for the effect of the thermal
motion on both the integrated cross-section as well as secondary angle and
energy distributions. For integrated cross-sections, it is possible to calculate
thermally-averaged cross-sections by applying a kernel Doppler broadening
algorithm to data at 0 K (or some temperature lower than the desired
temperature). The most ubiquitous algorithm for this purpose is the [SIGMA1]_
method developed by Red Cullen and subsequently refined by others. This method
is used in the NJOY_ and PREPRO_ data processing codes.

The effect of thermal motion on secondary angle and energy distributions can be
accounted for on-the-fly in a Monte Carlo simulation. We must first qualify
where it is actually used however. All threshold reactions are treated as being
independent of temperature, and therefore they are not Doppler broadened in NJOY
and no special procedure is used to adjust the secondary angle and energy
distributions. The only non-threhold reactions with secondary neutrons are
elastic scattering and fission. For fission, it is assumed that neutrons are
emitted isotropically (this is not strictly true, but is nevertheless a good
approximation). This leaves only elastic scattering that needs a special thermal
treatment for secondary distributions.

Fortunately, it is possible to directly sample the velocity of the target
nuclide and then use it directly in the kinematic calculations. However, this
calculation is a bit more nuanced than it might seem at first glance. One might
be tempted to simply sample a Maxwellian distribution for the velocity of the
target nuclide.  Careful inspection of equation :eq:`freegas1` however tells us
that target velocities that produce relative velocities which correspond to high
cross sections will have a greater contribution to the effection reaction
rate. This is most important when the velocity of the incoming neutron is close
to a resonance. For example, if the neutron's velocity corresponds to a trough
in a resonance elastic scattering cross-section, a very small target velocity
can cause the relative velocity to correpond to the peak of the resonance, thus
making a disproportionate contribution to the reaction rate. The conclusion is
that if we are to sample a target velocity in the Monte Carlo code, it must be
done in such a way that preserves the thermally-averaged reaction rate as per
equation :eq:`freegas`.

The method by which most Monte Carlo codes sample the target velocity for use in
elastic scattering kinematics is outlined in detail by [Gelbard]_. The
derivation here largely follows that of Gelbard. The first assumption we can
make is that the velocity distribution for the thermal motion is isotropic, i.e.

.. math::
    :label: freegas2

    p(\mathbf{v}_T) d\mathbf{v}_T = \frac{1}{4\pi} p(v_T) dv_T d\mu d\phi

With this assumption, we can now rewrite equation :eq:`freegas1` as

.. math::
    :label: freegas3

    v_n \sigma (v_n, T) = \frac{1}{2} \int_{-1}^1 d\mu \int\limits_{v_r > 0}
    dv_T v_r \sigma (v_r, 0) p(v_T)

after integrating over :math:`d\phi`. To change the outer variable of
integration from :math:`\mu` to :math:`v_r`, we can establish a relation between
these variables based on the law of cosines.

.. math::
    :label: lawcosine

    2 v_n v_T \mu = v_n^2 + v_T^2 - v_r^2

The probability distribution for the magnitude of the velocity of the target
nucleus and the angle between the neutron and target velocity is

.. math::
    :label: freegas4

    P(v_T, \mu) = \frac{\sigma (v_r, 0) v_r P(v_T)}{2 \sigma (v_n, T) v_n}

It is normally assumed that :math:`\sigma (v_r, 0)` is constant over the range
of relative velocities of interest. This is a good assumption for almost all
cases since the elastic scattering cross section varies slowly with velocity for
light nuclei, and for heavy nuclei where large variations can occur due to
resonance scattering, the moderating effect is rather small. Nonetheless, this
assumption can cause incorrect answers in systems with U-238 where the low-lying
resonances can cause a significant amount of upscatter that would be ignored by
this assumption.

With this (sometimes incorrect) assumption, we see that the probability
distribution is proportional to

.. math::
    :label: freegas5

    P(v_T, \mu) \propto v_r P(v_T) = | v_n - v_T | P(v_T)

We can divide this probability distribution into two parts as such:

.. math::
    :label: freegas6

    P(v_T, \mu) &= f_1(v_T, \mu) f_2(v_T) \\
    f_1(v_T, \mu) &= \frac{| v_n - v_T |}{C (v_n + v_T)} \\
    f_2(v_T) &= (v_n + v_T) P(v_T)

where :math:`C = \int dv_T \sigma v_r P(v_T)`. In general, any probability
distribution function of the form :math:`p(x) = f_1(x) f_2(x)` with
:math:`f_1(x)` bounded can be sampled by sampling :math:`x_s` from the
distribution

.. math::
    :label: freegas7

    \frac{f_2(x)}{\int f_2(x) dx}

and accepting it with probability

.. math::
    :label: freegas8

    \frac{f_1(x_s)}{\max f_1(x)}

It is normally assumed that the velocity distribution of the target nucleus
assumes a Maxwellian distribution in velocity.

------------
|sab| Tables
------------

For neutrons with thermal energies, generally less than 4 eV, the kinematics of
scattering can be affected by chemical binding and crystalline effects of the
target molecule. If these effects are not accounted for in a simulation, the
reported results may be highly inaccurate. There is no general analytic
treatment for the scattering kinematics at low energies, and thus when nuclear
data is processed for use in a Monte Carlo code, special tables are created that
give altered cross-sections and secondary angle/energy distributions for thermal
scattering. These tables are mainly used for moderating materials such as light
or heavy water, graphite, hydrogen in ZrH, beryllium, etc.

The theory behind |sab| is rooted in quantum mechanics and is quite
complex. Those interested in first principles derivations for formulae relating
to |sab| tables should be referred to the excellent books by [Williams]_ and
[Squires]_. For our purposes here, we will focus only on the use of already
processed data as it appears in the ACE format.

Each |sab| table can contain the following:

- Thermal inelastic scattering cross section
- Thermal elastic scattering cross section
- Correlated energy-angle distributions for thermal inelastic and elastic
  scattering

Note that when we refer to "inelastic" and "elastic" scattering now, we are
actually using these terms with respect to the *scattering system*. Thermal
inelastic scattering means that the scattering system is left in an excited
state, not any particular nucleus as is the case in inelastic level
scattering. In a crystalline material, the excitation could be the production of
phonons. In a molecule, it could be the excitation of rotational or vibrational
modes.

Both thermal elastic and thermal inelastic scattering are generally divided into
incoherent and coherent parts. Coherent elastic scattering refers to scattering
in crystalline solids like graphite or beryllium. These cross-sections are
characterized by the presence of "Bragg edges" that relate to the crystal
structure of the scattering material. Incoherent elastic scattering refers to
scattering in hydrogenous solids such as polyethylene. As it occurs in ACE data,
thermal inelastic scattering includes both coherent and incoherent effects and
is dominant for most other materials including hydrogen in water.

Calculating Integrated Cross Sections
-------------------------------------

The first aspect of using |sab| tables is calculating cross-sections to replace
the data that would normally appear on the incident neutron data, which do not
account for thermal binding effects. For incoherent elastic and inelastic
scattering, the cross-sections are stored as linearly interpolable functions on
a specified energy grid. For coherent elastic data, the cross section can be
expressed as

.. math::
    :label: coherent-elastic-xs

    \sigma(E) = \frac{\sigma_c}{E} \sum_{E_i < E} f_i e^{-4WE_i}.

where :math:`\sigma_c` is the effective bound coherent scattering cross section,
:math:`W` is the effective Debye-Waller coefficient, :math:`E_i` are the
energies of the Bragg edges, and :math:`f_i` are related to crystallographic
structure factors. Since the functional form of the cross-section is just 1/E
and the proportionality constant changes only at Bragg edges, the
proportionality constants are stored and then the cross-section can be
calculated analytically based on equation :eq:`coherent-elastic-xs`.

Outgoing Angle for Coherent Elastic Scattering
----------------------------------------------

The other aspect of using |sab| tables is determining the outgoing energy and
angle of the neutron after scattering. For incoherent and coherent elastic
scattering, the energy of the neutron does not actually change, but the angle
does change. For coherent elastic scattering, the angle will depend on which
Bragg edge scattered the neutron. The probability that edge :math:`i` will
scatter then neutron is given by

.. math::
    :label: coherent-elastic-probability

    \frac{f_i e^{-4WE_i}}{\sum_j f_j e^{-4WE_j}}.

After a Bragg edge has been sampled, the cosine of the angle of scattering is
given analytically by

.. math::
    :label: coherent-elastic-angle

    \mu = 1 - \frac{E_i}{E}

where :math:`E_i` is the energy of the Bragg edge that scattered the neutron. 

Outgoing Angle for Incoherent Elastic Scattering
------------------------------------------------

For incoherent elastic scattering, the probability distribution for the cosine
of the angle of scattering is represent as a series of equally-likely discrete
cosines :math:`\mu_{i,j}` for each incoming energy :math:`E_i` on the thermal
elastic energy grid. First the outgoing angle bin :math:`j` is sampled. Then, if
the incoming energy of the neutron satisfies :math:`E_i < E < E_{i+1}` the final
cosine is

.. math::
    :label: incoherent-elastic-angle

    \mu = \mu_{i,j} + f (\mu_{i+1,j} - \mu_{i,j})

where the interpolation factor is defined as

.. math::
    :label: sab-interpolation-factor

    f = \frac{E - E_i}{E_{i+1} - E_i}.

Outgoing Energy and Angle for Inelastic Scattering
--------------------------------------------------

On each |sab| table, there is a correlated angle-energy secondary distribution
for neutron thermal inelastic scattering. While the documentation for the ACE
format implies that there are a series of equiprobably outgoing energies, the
outgoing energies may have non-uniform probability distribution. In particular,
if the thermal data were processed with :math:`iwt = 0` in NJOY, then the first
and last outgoing energies have a relative probability of 1, the second and
second to last energies have a relative probability of 4, and all other energies
have a relative probability of 10. The procedure to determine the outgoing
energy and angle is as such. First, the inteprolation factor is determined from
equation :eq:`sab-interpolation-factor`. Then, an outgoing energy bin is sampled
either from a uniform distribution or from a skewed distribution as
discussed. The outgoing energy is then interpolated between values corresponding
to neighboring incoming energies:

.. math::
    :label: inelastic-energy

    E = E_{i,j} + f (E_{i+1,j} - E_{i,j})

where :math:`E_{i,j}` is the j-th outgoing energy corresponding to the i-th
incoming energy. For each combination of incoming and outgoing energies, there
is a series equiprobable outgoing cosines. An outgoing cosine bin is sampled
uniformly and then the final cosine is interpolated on the incoming energy grid:

.. math::
    :label: inelastic-angle

    \mu = \mu_{i,j,k} + f (\mu_{i+1,j,k} - \mu_{i,j,k})

where :math:`\mu_{i,j,k}` is the k-th outgoing cosine corresponding to the j-th
outgoing energy and the i-th incoming energy. 

----------------------------------------------
Unresolved Resonance Region Probability Tables
----------------------------------------------

In the unresolved resonance energy range, resonances may be so closely spaced
that it is not possible for experimental measurements to resolve all
resonances. To properly account for self-shielding in this energy range, OpenMC
uses the probability table method [Levitt]_. For most thermal reactors, the use
of probability tables will not significantly affect problem results. However,
for some fast reactors and other problems with an appreciable flux spectrum in
the unresolved resonance range, not using probability tables may lead to
incorrect results.

Probability tables in the ACE format are generated from the UNRESR module in
NJOY following the method of Levitt. A similar method employed for the RACER and
MC21_ Monte Carlo codes is described in [Sutton]_. For the discussion here, we
will focus only on use of the probability table table as it appears in the ACE
format.

Each probability table for a nuclide contains the following information at a
number of incoming energies within the unresolved resonance range:

- Cumulative probabilities for cross section bands
- Total cross section (or factor) in each band
- Elastic scattering cross section (or factor) in each band
- Fission cross section (or factor) in each band
- :math:`(n,\gamma)` cross section (or factor) in each band
- Neutron heating number (or factor) in each band

It should be noted that unresolved resonance probability tables affect only
integrated cross sections and no extra data need be given for secondary
angle/energy distributions. Secondary distributions for elastic and inelastic
scattering would be specified whether or not probability tables were present.

The procedure for determining cross sections in the unresolved range using
probability tables is as follows. First, the bounding incoming energies are
determined, i.e. find :math:`i` such that :math:`E_i < E < E_{i+1}`. We then
sample a cross section band :math:`j` using the cumulative probabilities for
table :math:`i`. This allows us to then calculate the elastic, fission, and
capture cross-sections from the probability tables interpolating between
neighboring incoming energies. If interpolation is specified, then
the cross sections are calculated as

.. math::
    :label: ptables-linlin

    \sigma = \sigma_{i,j} + f (\sigma_{i+1,j} - \sigma{i,j})

where :math:`f` is the interpolation factor defined in the same manner as
:eq:`sab-interpolation-factor`. If logarithmic interpolation is specified, the
cross sections are calculated as

.. math::
    :label: ptables-loglog

    \sigma = \exp \left ( \log \sigma_{i,j} + f \log
    \frac{\sigma_{i+1,j}}{\sigma_{i,j}} \right )

where the interpolation factor is now defined as

.. math::
    :label: log-interpolation-factor

    f = \frac{\log \frac{E}{E_i}}{\log \frac{E_{i+1}}{E_i}}

A flag is also present in the probability table that specifies whether an
inelastic cross section should be calculated. If so, this is done from a normal
reaction cross section (either MT=51 or a special MT). Finally, if the
cross-sections defined are above are specified to be factors and not true
cross-sections, they are multiplied by the underlying smooth cross section in
the unresolved range to get the actual cross sections. Lastly, the total cross
section is calculated as the sum of the elastic, fission, capture, and inelastic
cross sections.

.. _NJOY: http://t2.lanl.gov/codes.shtml

.. _PREPRO: http://www-nds.iaea.org/ndspub/endf/prepro/

.. _MC21: http://www.osti.gov/bridge/servlets/purl/903083-HT5p1o/903083.pdf

.. [SIGMA1] Dermett E. Cullen and Charles R. Weisbin, "Exact Doppler Broadening
   of Tabulated Cross Sections," *Nucl. Sci. Eng.*, **60**, pp. 199-229 (1976).

.. [Gelbard] Ely M. Gelbard, "Epithermal Scattering in VIM," FRA-TM-123, Argonne
   National Laboratory (1979).

.. [Williams] M. M. R. Williams, *The Slowing Down and Thermalization of
   Neutrons*, North-Holland Publishing Co., Amsterdam (1966).

.. [Squires] G. L. Squires, *Introduction to the Theory of Thermal Neutron
   Scattering*, Cambridge University Press (1978).

.. [Levitt] Leo B. Levitt, "The Probability Table Method for Treating Unresolved
   Neutron Resonances in Monte Carlo Calculations," *Nucl. Sci. Eng.*, **49**,
   pp. 450-457 (1972).

.. [Sutton] Thomas M. Sutton and Forrest B. Brown, "Implementation of
   the Probability Table Method in a Continuous-Energy Monte Carlo Code System,"
   *Proc. International Conf. on the Physics of Nucl. Sci. and Technology*,
   October 5-8, Long Island, New York (1998).

.. |sab| replace:: S(:math:`\alpha,\beta`)
