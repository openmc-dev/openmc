"""Persistent, process-safe history reservations around the existing Runner.

Cached observations consume no new reservation. Failed/interrupted attempts
retain their full reservation, conservatively bounding transport expenditure.
No transport or nuclear-data logic is reimplemented here.
"""
from dataclasses import asdict
from pathlib import Path
import json
import sqlite3
import h5py
import numpy as np
from sandbox.analysis_codes.runner import Runner
from sandbox.analysis_codes.storage import sample_id


class HistoryBudget:
    def __init__(self, directory, cap):
        self.path=Path(directory)/'history_budget.sqlite'
        self.path.parent.mkdir(parents=True,exist_ok=True)
        self.cap=int(cap)
        with sqlite3.connect(self.path,timeout=60) as db:
            db.execute('BEGIN IMMEDIATE')
            db.execute('CREATE TABLE IF NOT EXISTS metadata (initialized INTEGER)')
            db.execute('CREATE TABLE IF NOT EXISTS attempts (sample TEXT, histories INTEGER)')
            if db.execute('SELECT count(*) FROM metadata').fetchone()[0]==0:
                # Import pre-ledger checkpoints once. New workers see initialized
                # metadata and never scan the entire dataset independently.
                for file in (Path(directory)/'evaluations').glob('*.h5'):
                    with h5py.File(file,'r') as saved:
                        record=json.loads(saved.attrs['record'])
                    db.execute('INSERT INTO attempts VALUES (?,?)',
                        (record['identity'],record['particles']*record['batches']))
                db.execute('INSERT INTO metadata VALUES (1)')
            db.commit()
        if self.used()>self.cap:
            raise ValueError('Already-completed/reserved histories exceed the requested cap')

    def used(self):
        with sqlite3.connect(self.path,timeout=60) as db:
            return int(db.execute('SELECT coalesce(sum(histories),0) FROM attempts').fetchone()[0])

    def reserve(self, identity, histories):
        with sqlite3.connect(self.path,timeout=60) as db:
            db.execute('BEGIN IMMEDIATE')
            used=int(db.execute('SELECT coalesce(sum(histories),0) FROM attempts').fetchone()[0])
            if used+histories>self.cap:
                raise RuntimeError(f'History budget exhausted: {used}+{histories}>{self.cap}')
            db.execute('INSERT INTO attempts VALUES (?,?)',(identity,int(histories)))
            db.commit()


class BudgetedRunner(Runner):
    """Use existing physics/extraction/rollback with atomic history reservations."""
    def __init__(self,spec,output,cap):
        super().__init__(spec,output)
        self.history_budget=HistoryBudget(output,cap)

    def evaluate(self,values,*,stage,seed,collect_batches=False):
        # Reproduce the base Runner's identity solely to distinguish cache hits.
        payload=dict(provenance=self.provenance,
            parameters=[asdict(p) for p in self.parameters],
            values=np.asarray(values).tolist(),stage=stage,seed=seed,
            particles=self.spec.model.settings.particles,
            batches=self.spec.model.settings.batches)
        if collect_batches:payload['collect_batches']=True
        identity=sample_id(payload)
        existing=self.dataset.load(identity)
        if existing is None or existing['status']!='success':
            self.history_budget.reserve(identity,payload['particles']*payload['batches'])
        return super().evaluate(values,stage=stage,seed=seed,collect_batches=collect_batches)
