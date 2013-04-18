.. _usersguide_beginners:

============================
A Beginner's Guide to OpenMC
============================

--------------------
What does OpenMC do?
--------------------

In a nutshell, OpenMC simulates neutrons moving around randomly in a `nuclear
reactor`_ (or other fissile system). This is what's known as `Monte Carlo`_
simulation. Neutrons are important in nuclear reactors because they are the
particles that induce `fission`_ in uranium and other nuclides. Knowing the
behavior of neutrons allows you to figure out how often and where fission
occurs. The amount of energy released is then directly proportional to the
fission reaction rate since most heat is produced by fission. By simulating many
neutrons (millions or billions), it is possible to determine the average
behavior of these neutrons (or the behavior of the energy produced or any other
quantity one is interested in) very accurately.

Using Monte Carlo methods to determine the average behavior of various physical
quantities in a nuclear reactor is quite different from other means of solving
the same problem. The other class of methods for determining the behavior of
neutrons and reactions rates in a reactor is so-called `deterministic`_
methods. In these methods, the starting point is not randomly simulating
particles but rather writing an equation that describes the average behavior of
the particles. The equation that describes the average behavior of neutrons is
called the `neutron transport`_ equation. This equation is a seven-dimensional
equation (three for space, three for velocity, and one for time) and is very
difficult to solve directly. For all but the simplest problems, it is necessary
to make some sort of `discretization`_. As an example, we can divide up all
space into small sections which are homogeneous and then solve the equation on
those small sections. After these discretizations and various approximations,
one can arrive at forms that are suitable for solution on a computer. Among
these are discrete ordinates, method of characteristics, finite-difference
diffusion, and nodal methods.

So why choose Monte Carlo over deterministic methods? Each method has its pros
and cons. Let us first take a look at few of the salient pros and cons of
deterministic methods:

- **Pro**: Depending on what method is used, solution can be determined very
  quickly.

- **Pro**: The solution is a global solution, i.e. we know the average behavior
  everywhere.

- **Pro**: Once the problem is converged, the solution is known.

- **Con**: If the model is complex, it is necessary to do sophisticated mesh
  generation.

- **Con**: It is necessary to generate multi-group cross sections which requires
  knowing the solution *a priori*.

Now let's look at the pros and cons of Monte Carlo methods:

- **Pro**: No mesh generation is required to build geometry. By using
  `constructive solid geometry`_, it's possible to build arbitrarily complex
  models with curved surfaces.

- **Pro**: Monte Carlo methods can be used with either continuous-energy or
  multi-group cross sections.

- **Pro**: Running simulations in parallel is conceptually very simple.

- **Con**: Because they related on repeated random sampling, they are
  computationally very expensive.

- **Con**: A simulation doesn't automatically give you the global solution
  everywhere -- you have to specifically ask for those quantities you want.

- **Con**: Even after the problem is converged, it is necessary to simulate
  many particles to reduce stochastic uncertainty.

Because fewer approximations are made in solving a problem by the Monte Carlo
method, it is often seen as a "gold standard" which can be used as a benchmark
for a solution of the same problem by deterministic means. However, it comes at
the expense of a potentially longer simulation.

-----------------
How does it work?
-----------------

In order to do anything, the code first needs to have a model of some problem of
interest. This could be a nuclear reactor or any other physical system with
fissioning material. You, as the code user, will need to describe the model so
that the code can do something with it. A basic model consists of a few things:

- A description of the geometry -- the problem should be split up into regions
  of homogeneous material.
- For each different material in the problem, a description of what nuclides are
  in the material and at what density.
- Various parameters telling the code how many particles to simulate and what
  options to use.
- A list of different physical quantities that the code should return at the end
  of the simulation. Remember, in a Monte Carlo simulation, if you don't ask for
  anything, it will not give you any answers (other than a few default
  quantities).

-----------------------
What do I need to know?
-----------------------

If you are starting to work with OpenMC, there are a few things you should be
familiar with. Whether you plan on working in Linux, Mac OS X, or Windows, you
should be comfortable working in a command line environment. There are many
resources online for learning command line environments. If you are using Linux
or Mac OS X (also Unix-derived), `this tutorial
<http://www.ee.surrey.ac.uk/Teaching/Unix/>`_ will help you get acquainted with
commonly-used commands. It is also helpful to be familiar with `Python
<http://www.python.org/>`_, as most of the post-processing utilities provided
with OpenMC rely on it for data manipulation and results visualization.

OpenMC uses a version control software called `git`_ to keep track of changes to
the code, document bugs and issues, and other development tasks. While you don't
necessarily have to have git installed in order to download and run OpenMC, it
makes it much easier to receive updates if you do have it installed and have a
basic understanding of how it works. There are a list of good `git tutorials`_
at the git documentation website. The `OpenMC source code`_ and documentation
are hosted at `GitHub`_. In order to receive updates to the code directly,
submit `bug reports`_, and perform other development tasks, you may want to sign
up for a free account on GitHub. Once you have an account, you can follow `these
instructions <http://help.github.com/set-up-git-redirect>`_ on how to set up
your computer for using GitHub.

If you are new to nuclear engineering, you may want to review the NRC's `Reactor
Concepts Manual`_. This manual describes the basics of nuclear power for
electricity generation, the fission process, and the overall systems in a
pressurized or boiling water reactor. Another resource that is a bit more
technical than the Reactor Concepts Manual but still at an elementary level is
the DOE Fundamentals Handbook on Nuclear Physics and Reactor Theory `Volume I`_
and `Volume II`_. You may also find it helpful to review the following terms:

- `Neutron cross section`_
- `Effective multiplication factor`_
- `Flux`_

.. _nuclear reactor: http://en.wikipedia.org/wiki/Nuclear_reactor
.. _Monte Carlo: http://en.wikipedia.org/wiki/Monte_Carlo_method
.. _fission: http://en.wikipedia.org/wiki/Nuclear_fission
.. _deterministic: http://en.wikipedia.org/wiki/Deterministic_algorithm
.. _neutron transport: http://en.wikipedia.org/wiki/Neutron_transport
.. _discretization: http://en.wikipedia.org/wiki/Discretization
.. _constructive solid geometry: http://en.wikipedia.org/wiki/Constructive_solid_geometry
.. _git: http://git-scm.com/
.. _git tutorials: http://git-scm.com/documentation
.. _Reactor Concepts Manual: http://web.mit.edu/romano7/www/reactor_concepts.pdf
.. _Volume I: http://www.hss.doe.gov/nuclearsafety/techstds/docs/handbook/h1019v1.pdf
.. _Volume II: http://www.hss.doe.gov/nuclearsafety/techstds/docs/handbook/h1019v2.pdf
.. _OpenMC source code: https://github.com/mit-crpg/openmc
.. _GitHub: https://github.com/
.. _bug reports: https://github.com/mit-crpg/openmc/issues
.. _Neutron cross section: http://en.wikipedia.org/wiki/Neutron_cross_section
.. _Effective multiplication factor: http://en.wikipedia.org/wiki/Effective_multiplication_factor
.. _Flux: http://en.wikipedia.org/wiki/Neutron_flux

