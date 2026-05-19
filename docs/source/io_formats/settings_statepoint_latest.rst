.. _usersguide_statepoint_latest:

Running Statepoints
===================

OpenMC can be configured to keep running statepoint files that capture the state
of the most recently completed batches. This is useful when you want access to
results from the most recent batches even if a simulation is interrupted
unexpectedly, allowing easy restart or analysis.

To enable this feature, specify a negative batch number in the ``<batches>``
element within ``<state_point>``. A value of ``-N`` means "keep the last N
completed batches as running statepoint files."

Examples

Keep the last batch completed:

.. code-block:: xml

  <state_point>
    <batches>-1</batches>
  </state_point>

Keep the last two completed batches:

.. code-block:: xml

  <state_point>
    <batches>-2</batches>
  </state_point>

You can also mix positive and negative batches. Positive batches specify
explicit intervals to write statepoints, while negative values keep a rolling
window of recent batches:

.. code-block:: xml

  <state_point>
    <batches>10 20 30 -3</batches>
  </state_point>

In this case, statepoints are written at batches 10, 20, and 30, and the last 3
completed batches are also retained as running statepoints.

File Naming
-----------

Running statepoint files are named ``statepoint.running.<batch>.h5``, where
``<batch>`` is the batch number. For example:

- Batch 1: ``statepoint.running.1.h5``
- Batch 2: ``statepoint.running.2.h5``
- Batch 5: ``statepoint.running.5.h5``

Pruning
-------

When you specify ``-N``, OpenMC automatically keeps only the most recent ``N``
running statepoints. Older running statepoint files are automatically deleted
when a new batch completes and the number of running statepoints exceeds ``N``.

For example, with ``<batches>-2</batches>`` and if batches 1 through 5 complete:

- After batch 1 completes: ``statepoint.running.1.h5`` exists
- After batch 2 completes: ``statepoint.running.1.h5``, ``statepoint.running.2.h5`` exist
- After batch 3 completes: ``statepoint.running.2.h5``, ``statepoint.running.3.h5`` exist (batch 1 file deleted)
- After batch 4 completes: ``statepoint.running.3.h5``, ``statepoint.running.4.h5`` exist (batch 2 file deleted)
- After batch 5 completes: ``statepoint.running.4.h5``, ``statepoint.running.5.h5`` exist (batch 3 file deleted)

Remarks
-------

Running statepoints are written in addition to any permanent statepoint files
specified by positive batch numbers. Each running statepoint file is named with
its batch number, making it easy to identify which batch's results are
contained in each file.
