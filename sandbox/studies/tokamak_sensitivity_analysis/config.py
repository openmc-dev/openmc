"""Sensitivity priors and transport budgets for the DEMO-inspired tokamak.

Geometry and material settings are inherited from sandbox.tokomak_analysis/config.py.
Edit TOKAMAK_CONFIG below to change that model; this file controls the study.
"""
from dataclasses import asdict
from sandbox.tokomak_analysis.config import Config as TokamakConfig

STUDY_NAME = 'DEMO-inspired Li2O tokamak'
RESULTS_DIRECTORY = 'results_six_nuclides'
TOKAMAK_CONFIG = asdict(TokamakConfig(
    fidelity=2, particles=200, batches=10, threads=1,
    enable_mesh=False, enable_heating=False, enable_exposure=False,
))
PRODUCTION_CHANNELS = {
    'Li6': [105], 'Li7': list(range(52, 83)), 'Be9': [700, 701],
}
RESPONSE_LABEL = 'Verified Li6 + Li7 + Be9 tritium production [tritons/source neutron]'
TARGET_NUCLIDES = ('Li6', 'Li7', 'Be9', 'O16', 'Fe56')
PERTURBATION_LIMIT = 0.20
MORRIS_TRAJECTORIES = 10
MORRIS_LEVELS = 4
TOP_PER_FAMILY = 8
# Optional physics checks to retain in the PCE even when they do not rank in
# the per-family Morris selection. Each tuple is (nuclide, MT, family, mode).
RESERVED_PARAMETERS = ()
SCREEN_PARTICLES = 100
SCREEN_BATCHES = 10
HIGH_PARTICLES = 1000
HIGH_BATCHES = 10
TRAINING_SAMPLES = 128
VALIDATION_SAMPLES = 32
VERIFICATION_PARTICLES = 200
VERIFICATION_BATCHES = 10
SEED = 20261106
ANALYSIS_SEED = 20261107
WORKERS = 1
HIGH_WORKERS = 1
TRANSPORT_THREADS = 1
PAIRED_MORRIS = True
LOW_HISTORY_CAP = 500_000_000
HIGH_HISTORY_CAP = 100_000_000
VERIFICATION_HISTORY_CAP = 20_000_000
HISTORY_CAP = 650_000_000
PCE_BACKEND = 'numpy'
PCE_ORDERS = (2, 3)
PCE_CV_FOLDS = 4
PCE_HYPERBOLIC_Q = 0.75
PCE_MAX_INTERACTION = 2
PCE_RIDGE = 1.0e-6
PCE_SPARSE = True
PCE_KEEP_MAIN_EFFECTS = True
PCE_BASIS_LIMIT = 1200
PCE_DISTRIBUTION_SAMPLES = 100_000
PCE_MIN_R2 = 0.90
PCE_MAX_NOISE_MSE = 2.0
PCE_MIN_SIGNAL_TO_NOISE = 3.0
