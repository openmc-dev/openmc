"""One initialized native library per model, safe rollback and checkpointing."""

from contextlib import AbstractContextManager
from pathlib import Path
from tempfile import TemporaryDirectory
from dataclasses import asdict
import hashlib
import os
import subprocess
import time
import shutil
import numpy as np
import openmc
import openmc.lib as lib
from .perturbations import NativeState
from .storage import sample_id, Dataset


class Runner(AbstractContextManager):
    """Serialized native execution with full model-specific provenance."""

    def __init__(self, spec, output):
        self.spec = spec
        self.dataset = Dataset(Path(output) / "evaluations")
        self.poisoned = False

    def __enter__(self):
        self.temporary = TemporaryDirectory(prefix="openmc-sensitivity-")
        try:
            repository = Path(__file__).resolve().parents[2]
            self.spec.model.init_lib(
                directory=self.temporary.name,
                output=False,
                threads=self.spec.threads,
            )
            xml = b"".join(
                p.read_bytes()
                for p in sorted(Path(self.temporary.name).glob("*.xml"))
            )
            library = Path(
                os.environ.get(
                    "OPENMC_CROSS_SECTIONS",
                    openmc.config.get("cross_sections", ""),
                )
            ).resolve()
            self.provenance = dict(
                model=self.spec.name,
                model_sha256=hashlib.sha256(xml).hexdigest(),
                version=openmc.__version__,
                git_commit=subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], text=True, cwd=repository
                ).strip(),
                library=str(library),
                library_index_sha256=hashlib.sha256(
                    library.read_bytes()
                ).hexdigest(),
                threads=self.spec.threads,
                responses={
                    k: asdict(v) for k, v in self.spec.responses.items()
                },
            )
            library_binary = Path(lib._dll._name).resolve()
            self.provenance["compiled_library"] = str(library_binary)
            self.provenance["compiled_library_sha256"] = hashlib.sha256(
                library_binary.read_bytes()
            ).hexdigest()
            self.provenance["git_status"] = subprocess.check_output(
                ["git", "status", "--short"], text=True, cwd=repository
            )
            inputs = (
                self.dataset.directory.parent
                / "model_inputs"
                / self.provenance["model_sha256"]
            )
            inputs.mkdir(parents=True, exist_ok=True)
            for file in Path(self.temporary.name).glob("*.xml"):
                shutil.copy2(file, inputs / file.name)
            # Include source of the dirty runtime/framework, not only HEAD.
            files = [
                *Path(__file__).resolve().parent.glob("*.py"),
                repository / "src/secondary_data.cpp",
                repository / "src/nuclide.cpp",
                repository / "openmc/lib/perturbation.py",
            ]
            self.provenance["implementation_sha256"] = hashlib.sha256(
                b"".join(p.read_bytes() for p in sorted(files))
            ).hexdigest()
            return self
        except BaseException:
            self.spec.model.finalize_lib()
            self.temporary.cleanup()
            raise

    def bind(self, parameters):
        self.parameters = parameters
        self.state = NativeState(parameters)

    def __exit__(self, *args):
        try:
            self.spec.model.finalize_lib()
        finally:
            self.temporary.cleanup()

    def evaluate(
        self, values, *, stage: str, seed: int, collect_batches=False
    ) -> dict:
        """Run or resume an identical sample; failures never become responses."""
        if self.poisoned:
            raise RuntimeError("Runner unusable after failed restoration")
        payload = dict(
            provenance=self.provenance,
            parameters=[asdict(p) for p in self.parameters],
            values=np.asarray(values).tolist(),
            stage=stage,
            seed=seed,
            particles=self.spec.model.settings.particles,
            batches=self.spec.model.settings.batches,
        )
        if collect_batches:
            payload["collect_batches"] = True
        identity = sample_id(payload)
        cached = self.dataset.load(identity)
        if cached is not None and cached["status"] == "success":
            self.state.verify()
            return cached
        record = dict(
            **payload,
            identity=identity,
            status="failed",
            error="",
            responses={},
        )
        started = time.perf_counter()
        try:
            self.state.verify()
            lib.hard_reset()
            lib.settings.seed = seed
            candidate = self.state.candidates(values)
            with lib.nuclear_data_perturbation(**candidate):
                self.state.verify_candidates(candidate)
                weighted = any(r.score_weights is not None for r in self.spec.responses.values())
                if collect_batches or weighted:
                    names = [
                        name
                        for name in ("total", "radial")
                        if name in self.spec.responses
                    ]
                    if weighted:
                        names = list(self.spec.responses)
                    batches = {name: [] for name in names}
                    previous = {name: 0.0 for name in names}
                    previous_realizations = {name: 0 for name in names}
                    try:
                        for _ in range(lib.settings.get_batches()):
                            lib.next_batch()
                            for name in names:
                                tally = lib.tallies[
                                    self.spec.responses[name].tally_id
                                ]
                                realized = tally.num_realizations
                                if realized == previous_realizations[name]:
                                    continue  # Inactive batches are not response observations.
                                if realized != previous_realizations[name] + 1:
                                    raise RuntimeError("Unexpected tally realization increment")
                                previous_realizations[name] = realized
                                response = self.spec.responses[name]
                                accumulated = response.reduce_scores(
                                    tally.results[:, :, 1].copy(),
                                    len(tally.nuclides), len(tally.scores))
                                scale = self.spec.responses[name].scale(
                                    len(accumulated)
                                )
                                batches[name].append(
                                    (accumulated - previous[name]) * scale
                                )
                                previous[name] = accumulated
                    finally:
                        lib.simulation_finalize()
                    record["batch_responses"] = batches
                else:
                    lib.run()
                for name, response in self.spec.responses.items():
                    tally = lib.tallies[response.tally_id]
                    if response.score_weights is not None:
                        batch = np.asarray(record["batch_responses"][name], float)
                        if len(batch) < 2:
                            raise ValueError("Weighted production responses require >=2 batches")
                        # Contributions are summed WITHIN each batch. The
                        # resulting SE retains every within-history covariance.
                        record["responses"][name] = (
                            batch.mean(axis=0), batch.std(axis=0, ddof=1)/np.sqrt(len(batch)))
                    else:
                        mean = response.reduce_scores(tally.mean, len(tally.nuclides), len(tally.scores))
                        se = response.reduce_scores(tally.std_dev, len(tally.nuclides), len(tally.scores))
                        scale = response.scale(len(mean))
                        record["responses"][name] = (mean * scale, se * scale)
            self.state.verify()
            record["status"] = "success"
        except Exception as error:
            record["error"] = repr(error)
            try:
                self.state.verify()
            except Exception:
                self.poisoned = True
            self.dataset.save(identity, record)
            raise RuntimeError(
                f"Native sample failed and saved: {identity}: {error}"
            ) from error
        record["wall_seconds"] = time.perf_counter() - started
        self.dataset.save(identity, record)
        return record
