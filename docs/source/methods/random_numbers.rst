.. _methods_random_numbers:

========================
Random Number Generation
========================

In order to sample probability distributions, one must be able to produce random
numbers. The standard technique to do this is to generate numbers on the
interval :math:`[0,1)` from a deterministic sequence that has properties that
make it appear to be random, e.g. being uniformly distributed and not exhibiting
correlation between successive terms. Since the numbers produced this way are
not truly "random" in a strict sense, they are typically referred to as
pseudorandom numbers, and the techniques used to generate them are pseudorandom
number generators (PRNGs). Numbers sampled on the unit interval can then be
transformed for the purpose of sampling other continuous or discrete probability
distributions.

OpenMC currently uses PCG. A short description of LCG is included, since 
it is essential to understand PCG.

------------------------------
Linear Congruential Generators
------------------------------

There are a great number of algorithms for generating random numbers. One of the
simplest and commonly used algorithms is called a `linear congruential
generator`_. We start with a random number *seed* :math:`\xi_0` and a sequence
of random numbers can then be generated using the following recurrence relation:

.. math::
    :label: lcg

    \xi_{i+1} = g \xi_i + c \mod M

where :math:`g`, :math:`c`, and :math:`M` are constants. The choice of these
constants will have a profound effect on the quality and performance of the
generator, so they should not be chosen arbitrarily. As Donald Knuth stated in
his seminal work *The Art of Computer Programming*, "random numbers should not
be generated with a method chosen at random. Some theory should be used."
Typically, :math:`M` is chosen to be a power of two as this enables :math:`x
\mod M` to be performed using the bitwise AND operator with a bit mask. The
constants for the linear congruential generator used by default in OpenMC are
:math:`g = 2806196910506780709`, :math:`c = 1`, and :math:`M = 2^{63}` (see
`L'Ecuyer`_).

Skip-ahead Capability
---------------------

One of the important capabilities for a random number generator is to be able to
skip ahead in the sequence of random numbers. Without this capability, it would
be very difficult to maintain reproducibility in a parallel calculation. If we
want to skip ahead :math:`N` random numbers and :math:`N` is large, the cost of
sampling :math:`N` random numbers to get to that position may be prohibitively
expensive. Fortunately, algorithms have been developed that allow us to skip
ahead in :math:`O(\log_2 N)` operations instead of :math:`O(N)`. One algorithm
to do so is described in a paper by Brown_. This algorithm relies on the following
relationship:

.. math::
    :label: lcg-skipahead

    \xi_{i+k} = g^k \xi_i + c \frac{g^k - 1}{g - 1} \mod M

Note that equation :eq:`lcg-skipahead` has the same general form as equation :eq:`lcg`, so
the idea is to determine the new multiplicative and additive constants in
:math:`O(\log_2 N)` operations.


--------------------------------
Permuted Congruential Generators
--------------------------------

A permuted congruential generator (PCG) aims to improve statistical quality 
of a LCG by using permutation functions. The general algorithm consists of 
two steps:

1. advance the LCG generator according to :eq:`pcg_lcg`,
2. output the LCG state "scrambled" by permutation function.



OpenMC uses the PCG-RXS-M-XS variant with 64-bit state and 
64-bit output. The code for permutation functiom is taken 
from `PCG GitHub`_. The exact algorithm follows.

1. Generate LCG output :math:`\xi_{i+1}`

  .. math::
      :label: pcg_lcg
  
      \xi_{i+1} = g \xi_i + c

  where :math:`g=6364136223846793005` and :math:`c=1442695040888963407`.

2. Apply random xorshift. 
    * First, take the upper 5 bits of :math:`\xi_{i+1}` 
      using the `arithmetic right shift operator`_ by 59 positions.
    * Second, add :math:`5` and shift seed by that number of bits
    * Third, perform `bitwise XOR`_ with the upper bits and original :math:`\xi_{i+1}`.

#. multiply by :math:`6364136223846793005`,
#. apply xorshift - similarly to step 2, use `arithmetic right shift operator`_ 
   to take upper 21 bits and perform `bitwise XOR`_.

As this might be difficult to imagine, let's add an example. 

* Let's assume :math:`\xi_{i} = 1`.
* :math:`\xi_{i+1} = 7806831264735756412` according to :eq:`pcg_lcg`, which 
  is ``0110 1100 0101 0111 0110 1111 1010 1100 0100 0011 1111 1101 0000 0000 0111 1100`` in bit representation.
* After performing the bit shift by 59 positions the number in bits is 
  ``0 1101`` or :math:`13`, when represented as integer.
* Adding 5, we have to shift :math:`\xi_{i+1}` by :math:`18` bits, which yields
  ``01 1011 0001 0101 1101 1011 1110 1011 0001 0000 1111 1111`` or :math:`29780697878783` when represented as integer. 
* Perform `bitwise XOR`_ with shifted bits and original :math:`\xi_{i+1}`::

    0110 1100 0101 0111 0110 1111 1010 1100 0100 0011 1111 1101 0000 0000 0111 1100
    0000 0000 0000 0000 0001 1011 0001 0101 1101 1011 1110 1011 0001 0000 1111 1111
    -------------------------------------------------------------------------------
    0110 1100 0101 0111 0111 0100 1011 1001 1001 1000 0001 0110 0001 0000 1000 0011 

  The resulting number is :math:`7806836819539398787` as integer.
* After multiplication we get :math:`7806836819539398787 \cdot 6364136223846793005 = 13112265920887772427` 
  which is ``1011 0101 1111 1000 0010 0000 0110 0010 0001 1110 1111 1001 1101 0101 0000 1011`` as bits
* Shifted right by 43 positions is ``1 0110 1011 1111 0000 0100`` or :math:`1490692` as int.
* Finally, there is another XOR::

    1011 0101 1111 1000 0010 0000 0110 0010 0001 1110 1111 1001 1101 0101 0000 1011
    0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0001 0110 1011 1111 0000 0100
    -------------------------------------------------------------------------------
    1011 0101 1111 1000 0010 0000 0110 0010 0001 1110 1110 1111 0110 1010 0000 1111 

  And the output value as integer is :math:`13112265920887089679`.
* Convert the value to double from interval :math:`[0, 1)` as 
  :math:`13112265920887089679\cdot 2^{-64} = 0.710817`, which is the output of the generator.


For elaborated description, see `O'Neill`_.

**Advantages of PCG over LCG include**

* increased statistical quality - measured by statistical tests from BigCrush library,
* small performance burden compared to LCG.


.. only:: html

   .. rubric:: References


.. _L'Ecuyer: https://doi.org/10.1090/S0025-5718-99-00996-5
.. _Brown: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/anl-rn-arb-stride.pdf
.. _linear congruential generator: https://en.wikipedia.org/wiki/Linear_congruential_generator
.. _O'Neill: https://www.pcg-random.org/pdf/hmc-cs-2014-0905.pdf
.. _PCG GitHub: https://github.com/imneme/pcg-c/blob/83252d9c23df9c82ecb42210afed61a7b42402d7/include/pcg_variants.h#L188-L192
.. _arithmetic right shift operator: https://stackoverflow.com/a/141873/13224210
.. _bitwise XOR: https://www.learncpp.com/cpp-tutorial/bitwise-operators/
