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

All new features, enhancements, and bug fixes should be developed on a branch
that branches off of *develop*. When the feature is completed, a `pull request`_
is initiated on GitHub that is then reviewed by a committer. If the pull request
is satisfactory, it is then merged into *develop*. Note that a committer may not
review their own pull request (i.e., an independent code review is required).

Code Review Criteria
--------------------

In order to be considered suitable for inclusion in the *develop* branch, the
following criteria must be satisfied for all proposed changes:

- Changes have a clear purpose and are useful.
- Compiles and passes all tests under multiple build configurations (This is
  checked by Travis CI).
- If appropriate, test cases are added to regression or unit test suites.
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

1. Fork the main openmc repository from `openmc-dev/openmc`_. This will create a
   repository with the same name under your personal account. As such, you can
   commit to it as you please without disrupting other developers.

   .. image:: ../_images/fork.png

2. Clone your fork of OpenMC and create a branch that branches off of *develop*:

   .. code-block:: sh

       git clone --recurse-submodules git@github.com:yourusername/openmc.git
       cd openmc
       git checkout -b newbranch develop

3. Make your changes on the new branch that you intend to have included in
   *develop*. If you have made other changes that should not be merged back,
   ensure that those changes are made on a different branch.

4. Issue a pull request from GitHub and select the *develop* branch of
   openmc-dev/openmc as the target.

   At a minimum, you should describe what the changes you've made are and why
   you are making them. If the changes are related to an oustanding issue, make
   sure it is cross-referenced.

5. A committer will review your pull request based on the criteria
   above. Any issues with the pull request can be discussed directly on the pull
   request page itself.

6. After the pull request has been thoroughly vetted, it is merged back into the
   *develop* branch of openmc-dev/openmc.

Private Development
-------------------

While the process above depends on the fork of the OpenMC repository being
publicly available on GitHub, you may also wish to do development on a private
repository for research or commercial purposes. The proper way to do this is to
create a complete copy of the OpenMC repository (not a fork from GitHub). The
private repository can then either be stored just locally or in conjunction with
a private repository on Github (this requires a `paid plan`_). Alternatively,
`Bitbucket`_ offers private repositories for free. If you want to merge some
changes you've made in your private repository back to openmc-dev/openmc
repository, simply follow the steps above with an extra step of pulling a branch
from your private repository into a public fork.

.. _devguide_editable:

Working in "Development" Mode
-----------------------------

If you are making changes to the Python API during development, it is highly
suggested to install the Python API in development/editable mode using
pip_. From the root directory of the OpenMC repository, run:

.. code-block:: sh

    pip install -e .[test]

This installs the OpenMC Python package in `"editable" mode
<https://pip.pypa.io/en/stable/reference/pip_install/#editable-installs>`_ so
that 1) it can be imported from a Python interpreter and 2) any changes made are
immediately reflected in the installed version (that is, you don't need to keep
reinstalling it). While the same effect can be achieved using the
:envvar:`PYTHONPATH` environment variable, this is generally discouraged as it
can interfere with virtual environments.

.. _git: http://git-scm.com/
.. _GitHub: https://github.com/
.. _git flow: http://nvie.com/git-model
.. _valgrind: http://valgrind.org/
.. _style guide: https://docs.openmc.org/en/latest/devguide/styleguide.html
.. _pull request: https://help.github.com/articles/using-pull-requests
.. _openmc-dev/openmc: https://github.com/openmc-dev/openmc
.. _paid plan: https://github.com/plans
.. _Bitbucket: https://bitbucket.org
.. _pip: https://pip.pypa.io/en/stable/
