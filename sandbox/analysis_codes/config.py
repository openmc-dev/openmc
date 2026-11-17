"""Model registration and response semantics."""

from dataclasses import dataclass, field
from typing import Any
import numpy as np


@dataclass
class Response:
    """A tally response; optional volumes convert bin rates to rate densities.

    Without score_weights, select a single score and nuclide. Weighted sums
    require native per-batch scoring before estimating uncertainties. Mesh bins
    retain OpenMC's native filter order; regular mesh x varies fastest.
    """

    tally_id: int
    units: str = "reactions/source"
    volumes: Any = None
    coordinates: Any = None
    label: str = "reaction rate"
    score_weights: Any = None

    def reduce_scores(self, values, nuclides: int, scores: int) -> np.ndarray:
        """Reduce native nuclide/score bins before batch statistics.

        Weights have shape (nuclides, scores). This supports production-yield
        sums without dropping covariance. Native tally filter order is retained.
        """
        values = np.asarray(values, float).reshape(-1, nuclides * scores)
        if self.score_weights is None:
            if nuclides != 1 or scores != 1:
                raise ValueError("Response requires one score/nuclide or explicit score_weights")
            return values.ravel().copy()
        weights = np.asarray(self.score_weights, float)
        if (weights.shape != (nuclides, scores)
                or not np.all(np.isfinite(weights)) or np.any(weights < 0)):
            raise ValueError("Invalid nonnegative nuclide/score reduction weights")
        return values @ weights.ravel()

    def scale(self, size: int) -> np.ndarray:
        if self.volumes is None:
            return np.ones(size)
        volumes = np.broadcast_to(np.asarray(self.volumes, float), (size,))
        if not np.all(np.isfinite(volumes)) or np.any(volumes <= 0):
            raise ValueError("Volumes must be positive and finite")
        return 1 / volumes


@dataclass
class ModelSpec:
    """Named model with independently configured geometry, tallies and settings."""

    name: str
    model: Any
    responses: dict[str, Response] = field(default_factory=dict)
    target_nuclides: tuple[str, ...] = ("Li6",)
    active_ids: set[str] | None = None
    threads: int = 1

    def __post_init__(self):
        if not self.name or "/" in self.name or self.name in (".", ".."):
            raise ValueError("Model name must be a single path component")
        if not self.responses:
            for tally in self.model.tallies:
                self.responses[tally.name or str(tally.id)] = Response(
                    tally.id
                )
        if not self.responses:
            raise ValueError("Supply response tallies")
