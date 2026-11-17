"""Atomic per-evaluation checkpoints and compact HDF5 response arrays."""

import hashlib
import json
from pathlib import Path
import h5py
import numpy as np


def jsonable(value):
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, np.generic):
        return jsonable(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def save_json(path, value) -> None:
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(jsonable(value), indent=2, allow_nan=False)
    )
    temporary.replace(path)


def sample_id(payload: dict) -> str:
    """Hash full model/provenance/configuration/input identity."""
    return hashlib.sha256(
        json.dumps(jsonable(payload), sort_keys=True).encode()
    ).hexdigest()


class Dataset:
    """One independently atomic HDF5 file per run; metadata inside each file."""

    def __init__(self, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)

    def load(self, identity: str):
        path = self.directory / f"{identity}.h5"
        if not path.exists():
            return None
        with h5py.File(path, "r") as f:
            record = json.loads(f.attrs["record"])
            record["responses"] = {
                name: (f[name]["mean"][:], f[name]["se"][:]) for name in f
            }
        return record

    def save(self, identity: str, record: dict) -> None:
        path = self.directory / f"{identity}.h5"
        temporary = path.with_suffix(".tmp")
        metadata = {k: v for k, v in record.items() if k != "responses"}
        with h5py.File(temporary, "w") as f:
            f.attrs["record"] = json.dumps(jsonable(metadata))
            for name, (mean, se) in record.get("responses", {}).items():
                group = f.create_group(name)
                group.create_dataset("mean", data=mean, compression="gzip")
                group.create_dataset("se", data=se, compression="gzip")
        temporary.replace(path)
