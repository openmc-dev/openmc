.. _methods_physics:

=======
Physics
=======

------------------------------
Free Gas Scattering Kinematics
------------------------------

When a neutron scatters off of a nucleus, many times it is assumed that the
target nucleus is at rest. However, if the material is at a temperature greater
than 0 K, it will have motion associated with the thermal vibration. Thus, the
velocity of the neutrno relative to the target nucleus is in general not the
same as the velocity of the neutron entering the collision.

The affect of the thermal motion on the interaction probability can be written
as

.. math::
    :label: freegas1

    v_n \sigma (v_n, T) = \int_0^\infty d\mathbf{v}_T \sigma(v_r, 0)
    \mathbf{v}_r p(\mathbf{v}_T)
    
One assumption we can make here is that the velocity distribution for the
thermal motion is isotropic, i.e.

.. math::
    :label: freegas2

    p(\mathbf{v}_T) d\mathbf{v}_T = \frac{1}{4\pi} p(v_T) dv_T d\mu d\phi

With this assumption, we can now rewrite equation :eq:`freegas1` as

.. math::
    :label: freegas3

    v_n \sigma (v_n, T) = \frac{1}{2} \int_{-1}^1 d\mu \int\limits_{v_r > 0}
    v_r \sigma (v_r, 0) p(v_T) dv_T

To change the outer variable of integration from :math:`\mu` to :math:`v_r`, we
can establish a relation between these variables based on the law of cosines.

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
    f_1(v_T, \mu) &= \frac{| v_n - v_T |}{\hat{f_1} (v_n + v_T)} \\
    f_2(v_T) &= (v_n + v_T) P(v_T)

In general, any probability distribution function of the form :math:`p(x) =
f_1(x) f_2(x)` with :math:`f_1(x)` bounded can be sampled by sampling
:math:`x_s` from the distribution

.. math:: \frac{f_2(x)}{\int f_2(x) dx}

and accepting it with probability

.. math:: \frac{f_1(x_s)}{\max f_1(x)}

It is normally assumed that the velocity distribution of the target nucleus
assumes a Maxwellian distribution in velocity.
