.. _usersguide_statepoint_latest:

Statepoint Overwrite
====================

OpenMC can be configured to write a continuously-overwritten "running" statepoint
file after every finished batch. This is useful when you want the most recent
results available even if the job is interrupted unexpectedly.

To enable this behavior, set the ``overwrite_latest`` option inside the
``<state_point>`` element in the settings XML. The value may be either a
boolean or an integer:

- ``true``: keep a single file named ``statepoint.running.h5`` that is
  overwritten after every batch.
- ``false``: (default) do not write running statepoint files.
- ``N`` (integer &gt; 1): keep ``N`` rotating running statepoint files named
  ``statepoint.running.1.h5``, ``statepoint.running.2.h5``, ..., ``statepoint.running.N.h5``.

Examples

Keep a single running statepoint file:

.. code-block:: xml

  <state_point>
    <overwrite_latest>true</overwrite_latest>
  </state_point>

Keep the last two completed batches in rotating files:

.. code-block:: xml

  <state_point>
    <overwrite_latest>2</overwrite_latest>
  </state_point>

Remarks
-------

The running statepoints are written in addition to any permanent statepoint
files you requested via the ``<state_point><batches>...</batches></state_point>``
element. The rotating files form a circular buffer of the most recent
completed batches and are overwritten in batch order.
