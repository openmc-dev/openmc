.. _workflow:

====================
Development Workflow
====================

Anyone wishing to make contributions to OpenMC should be fully acquianted and
comfortable working with git_ and GitHub_. The primary means of modifying and
making contributions to OpenMC is through GitHub `pull requests`_. This is
what's known as a fork and pull development model. The steps for this are as
follows:

1. Fork the main openmc repository from `mit-crpg/openmc`_. This will create a
repository with the same name under your personal account. As such, you can
commit to it as you please without disrupting other developers.

2. Create a branch that you want merged back to `mit-crpg/openmc`_ and make
commits that you intend to go back. If you have made other changes that should
not be merged back, those changes should be on another branch.

3. Issue a pull request from GitHub and select the branch you want merged.

4. The OpenMC integration manager will review your pull request and make sure it
conforms to the :ref:`styleguide`, compiles correctly, runs in parallel
correctly, does not break other features in the code, etc. Any issues with the
pull request can be discussed directly on the pull request page itself.

5. After the pull request has been thoroughly vetted, it is merged back into
`mit-crpg/openmc`_.

While the process above depends on the fork of the OpenMC repository being
publicly available on GitHub, you may also wish to do development on a private
repository for research or commercial purposes. The proper way to do this is to
create a complete copy of the OpenMC repository (not a fork from GitHub). The
private repository can then either be stored just locally or in conjunction with
a private repository on Github (this requires a `paid plan`_). If you want to
merge some changes you've made in your private repository back to
`mit-crpg/openmc`_ repository, simply follow the steps above with an extra step
of pulling a branch from your private repository into your public fork.

.. _git: http://git-scm.com/
.. _GitHub: https://github.com/
.. _pull requests: https://help.github.com/articles/using-pull-requests
.. _mit-crpg/openmc: https://github.com/mit-crpg/openmc
.. _paid plan: https://github.com/plans
