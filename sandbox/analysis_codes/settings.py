"""Internal Morris settings; study config.py supplies every execution budget."""
from dataclasses import dataclass


@dataclass(frozen=True)
class Settings:
    trajectories: int = 20
    levels: int = 6
    particles: int = 10000
    batches: int = 10
    active: int = 8
    training: int = 256
    validation: int = 64
    seed: int = 20261020
    limit: float = .20
    orders: tuple = (2, 3)
    max_trajectories: int = 40
    trajectory_increment: int = 10
    max_histories: int | None = None
    workers: int = 8
    minimum: int = 20
    rank_correlation: float = .9
    top_overlap: float = .8
    max_relative_ci_width: float = .5
    bootstraps: int = 1000
    progress_bootstraps: int = 300
    convergence_bootstraps: int = 200
    confidence_level: float = .95
    noise_multiplier: float = 2.
    bootstrap_seed: int | None = None


FULL = Settings()