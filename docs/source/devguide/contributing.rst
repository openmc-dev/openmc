.. _devguide_contributing:

======================
Contributing to OpenMC
======================

Thank you for considering contributing to OpenMC! We look forward to welcoming
new members to the community and will do our best to help you get up to speed.
The purpose of this section is to document how the project is managed: how
contributions (bug fixes, enhancements, new features) are made, how they are
evaluated, who is permitted to merge pull requests, and what happens in the
event of disagreements. Once you have read through this section, the
:ref:`devguide_workflow` section outlines the actual mechanics of making a
contribution (forking, submitting a pull request, etc.).

The goal of our governance model is to:

- Encourage new contributions.
- Encourage contributors to remain involved.
- Avoid unnecessary processes and bureaucracy whenever possible.
- Create a transparent decision making process which makes it clear how
  contributors can be involved in decision making.

Overview
--------

OpenMC uses a liberal contribution model for project governance. Anyone involved
in development in a non-trivial capacity is given an opportunity to influence
the direction of the project. Project decisions are made through a
consensus-seeking process rather than by voting.

Terminology
-----------

- A *Contributor* is any individual creating or commenting on an issue or pull
  request.
- A *Committer* is a subset of contributors who are authorized to review and
  merge pull requests.
- The *TC* (Technical Committee) is a group of committers who have the authority
  to make decisions on behalf of the project team in order to resolve disputes.
- The *Project Lead* is a single individual who has the authority to make a final
  decision when the TC is unable to reach consensus.

Contribution Process
--------------------

Any change to the OpenMC repository must be made through a pull request (PR).
This applies to all changes to documentation, code, binary files, etc. Even long
term committers and TC members must use pull requests.

No pull request may be merged without being reviewed.

For non-trivial contributions, pull requests should not be merged for at least
36 hours to ensure that contributors in other timezones have time to review.
Consideration should also be given to weekends and other holiday periods to
ensure active committers have reasonable time to become involved in the
discussion and review process if they wish.

The default for each contribution is that it is accepted once no committer has
an objection. During review committers may also request that a specific
contributor who is most versed in a particular area review the PR before it can
be merged. Once all issues brought by committers are addressed it can be merged
by any committer.

In the case of an objection being raised in a pull request by another committer,
all involved committers should seek to arrive at a consensus by way of
addressing concerns being expressed through discussion, compromise on the
proposed change, or withdrawal of the proposed change.

If objections to a PR are made and committers cannot reach a consensus on how to
proceed, the decision is escalated to the TC. TC members should regularly
discuss pending contributions in order to find a resolution. It is expected that
only a small minority of issues be brought to the TC for resolution and that
discussion and compromise among committers be the default resolution mechanism.

Becoming a Committer
--------------------

All contributors who make a non-trivial contribution will be added as a
committer in a timely manner. Committers are expected to follow this policy.

TC Process
----------

The TC uses a "consensus seeking" process for issues that are escalated to the
TC. The group tries to find a resolution that has no objections among TC
members. If a consensus cannot be reached, the Project Lead has the ultimate
authority to make a final decision. It is expected that the majority of
decisions made by the TC are via a consensus seeking process and that the
Project Lead intercedes only as a last resort.

Resolution may involve returning the issue to committers with suggestions on how
to move forward towards a consensus.

Members can be added to the TC at any time. Any committer can nominate another
committer to the TC and the TC uses its standard consensus seeking process to
evaluate whether or not to add this new member. Members who do not participate
consistently at the level of a majority of the other members are expected to
resign.

In the event that the Project Lead resigns or otherwise steps down, the TC uses
a consensus seeking process to choose a new Project Lead.

Leadership Team
---------------

The TC consists of the following individuals:

- `Paul Romano <https://github.com/paulromano>`_
- `Sterling Harper <https://github.com/smharper>`_
- `Adam Nelson <https://github.com/nelsonag>`_
- `Benoit Forget <https://github.com/bforget>`_

The Project Lead is Paul Romano.