.. _devguide_workflow:

====================
Development Workflow
====================

Anyone wishing to make contributions to OpenMC should be fully acquianted and
comfortable working with git_ and GitHub_. We assume here that you have git
installed on your system, have a GitHub account, and have setup SSH keys to be
able to create/push to repositories on GitHub. 

Overview
--------

Development of OpenMC relies heavily on branching; specifically, we use a
branching model sometimes referred to as `git flow`_. If you plan to contribute
to OpenMC development, we highly recommend that you read the linked blog post to
get a sense of how the branching model works. There are two main branches that
always exist: *master* and *develop*. The *master* branch is a stable branch
that contains the latest release of the code. The *develop* branch is where any
ongoing development takes place prior to a release and is not guaranteed to be
stable. When the development team decides that a release should occur, the
*develop* branch is merged into *master*.

Trivial changes to the code may be committed directly to the *develop* branch by
a trusted developer. However, most new features should be developed on a branch
that branches off of *develop*. When the feature is completed, a `pull request`_
is initiated on GitHub that is then reviewed by a trusted developer. If the pull
request is satisfactory, it is then merged into *develop*. Note that a trusted
developer may not review their own pull request (i.e., an independent code
review is required).

Code Review Criteria
--------------------

In order to be considered suitable for inclusion in the *develop* branch, the
following criteria must be satisfied for all proposed changes:

- Changes have a clear purpose and are useful.
- Compiles under all conditions (MPI, OpenMP, HDF5, etc.).  This is checked as
  part of the test suite (see `test_compile.py`_).
- Passes the regression suite.
- If appropriate, test cases are added to regression suite.
- No memory leaks (checked with valgrind_).
- Conforms to the OpenMC `style guide`_.
- No degradation of performance or greatly increased memory usage. This is not a
  hard rule -- in certain circumstances, a performance loss might be acceptable
  if there are compelling reasons.
- New features/input are documented.
- No unnecessary external software dependencies are introduced.

Contributing
------------

Now that you understand the basic development workflow, let's discuss how an
individual to contribute to development. Note that this would apply to both new
features and bug fixes. The general steps for contributing are as follows:

1. Fork the main openmc repository from `mit-crpg/openmc`_. This will create a
   repository with the same name under your personal account. As such, you can
   commit to it as you please without disrupting other developers.

   .. image:: ../_images/fork.png

2. Clone your fork of OpenMC and create a branch that branches off of *develop*:

   .. code-block:: sh

       git clone git@github.com:yourusername/openmc.git
       cd openmc
       git checkout -b newbranch develop

3. Make your changes on the new branch that you intend to have included in
   *develop*. If you have made other changes that should not be merged back, 
   ensure that those changes are made on a different branch.

4. Issue a pull request from GitHub and select the *develop* branch of
   mit-crpg/openmc as the target.

   .. image:: ../_images/pullrequest.png

   At a minimum, you should describe what the changes you've made are and why
   you are making them. If the changes are related to an oustanding issue, make
   sure it is cross-referenced. A wise developer would also check whether their
   changes do indeed pass the regression test suite.

5. A trusted developer will review your pull request based on the criteria
   above. Any issues with the pull request can be discussed directly on the pull
   request page itself.

6. After the pull request has been thoroughly vetted, it is merged back into the
   *develop* branch of mit-crpg/openmc.

Private Development
-------------------

While the process above depends on the fork of the OpenMC repository being
publicly available on GitHub, you may also wish to do development on a private
repository for research or commercial purposes. The proper way to do this is to
create a complete copy of the OpenMC repository (not a fork from GitHub). The
private repository can then either be stored just locally or in conjunction with
a private repository on Github (this requires a `paid plan`_). Alternatively,
`Bitbucket`_ offers private repositories for free. If you want to merge some
changes you've made in your private repository back to mit-crpg/openmc
repository, simply follow the steps above with an extra step of pulling a branch
from your private repository into a public fork.

.. _git: http://git-scm.com/
.. _GitHub: https://github.com/
.. _git flow: http://nvie.com/git-model
.. _test_compile.py: https://github.com/mit-crpg/openmc/blob/develop/tests/test_compile/test_compile.py
.. _valgrind: http://valgrind.org/
.. _style guide: http://mit-crpg.github.io/openmc/devguide/styleguide.html
.. _pull request: https://help.github.com/articles/using-pull-requests
.. _mit-crpg/openmc: https://github.com/mit-crpg/openmc
.. _paid plan: https://github.com/plans
.. _Bitbucket: https://bitbucket.org
