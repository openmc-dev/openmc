"""Editable lithium benchmark and bounded study settings (units shown)."""
RADIUS_CM = 40.0
DENSITY_G_CM3 = 1.94
TEMPERATURE_K = 900.0
SOURCE_ENERGY_EV = 14.0e6
LI6_ATOM_FRACTION = None  # None: natural lithium; otherwise an isotope fraction.

TOP_PER_FAMILY = 8
MORRIS_TRAJECTORIES = 20
MORRIS_LEVELS = 6
SCREEN_PARTICLES = 10000
SCREEN_BATCHES = 10
HIGH_PARTICLES = 100000
HIGH_BATCHES = 10
TRAINING_SAMPLES = 512
VALIDATION_SAMPLES = 96
SEED = 20261103
HISTORY_CAP = 160_000_0000
LOW_HISTORY_CAP = 60_000_0000
HIGH_HISTORY_CAP = 80_000_0000
VERIFICATION_HISTORY_CAP = 20_000_0000
WORKERS = 8
ANALYSIS_SEED = 20261020  # CV folds and cheap PCE propagation; independent of transport seeds.

# Shared engine: run `python -m sandbox.analysis_codes FLiBe_analysis --run`.
STUDY_NAME = 'Li'
TARGET_NUCLIDES = None  # None: discover and screen every material nuclide.
TRANSPORT_THREADS = 1  # Per worker; WORKERS controls Morris process count.
RESULTS_DIRECTORY = 'results'  # Relative to this study directory.
PERTURBATION_LIMIT = 0.20  # Fractional bound: 0.10 = ±10%; finite 0 < limit < 1.
PAIRED_MORRIS = True  # Common seed within a trajectory; observed paired-batch SE.

# PCE backend: 'numpy' (weighted sparse/ridge) or 'openturns' (LARS/CLOO).
# OpenTURNS is optional: python -m pip install openturns
PCE_BACKEND = 'numpy'
PCE_ORDERS = (2, 3)
PCE_CV_FOLDS = 4
PCE_HYPERBOLIC_Q = 0.75
PCE_MAX_INTERACTION = 2
PCE_RIDGE = 1.e-6
PCE_SPARSE = True
PCE_KEEP_MAIN_EFFECTS = True
PCE_BASIS_LIMIT = 1200
PCE_DISTRIBUTION_SAMPLES = 100_000

# Verified ENDF/B-VIII.0 unit-triton-yield channels, locked to library hashes.
# Be9 MT105 is redundant: perturb/score its partial MT700 and MT701 instead.
PRODUCTION_CHANNELS = {'Li6': [105], 'Li7': list(range(52, 83))}
RESPONSE_LABEL = 'Total lithium tritium production [tritons/source neutron]'
VERIFICATION_PARTICLES = 10000
VERIFICATION_BATCHES = 20
# Reserve these native physics checks after top-n selection in EACH family.
RESERVED_PARAMETERS = [
    ('Li7', 2, 'xs', 0), ('Li7', 2, 'angle', 1), ('Li7', 2, 'angle', 2),
    ('Li6', 105, 'xs', 0), ('Li7', 56, 'xs', 0)
]

# Geometry/source controls. None places the source just outside the -z surface.
WORLD_MARGIN_CM = 1.0
SOURCE_POSITION_CM = None
SOURCE_STANDOFF_CM = 0.01
SOURCE_DIRECTION = (0.0, 0.0, 1.0)

# Morris uncertainty and convergence diagnostics (existing defaults).
MORRIS_MIN_TRAJECTORIES = 20
MORRIS_RANK_CORRELATION = 0.90
MORRIS_TOP_OVERLAP = 0.80
MORRIS_MAX_RELATIVE_CI_WIDTH = 0.50
MORRIS_BOOTSTRAPS = 1000
MORRIS_PROGRESS_BOOTSTRAPS = 300
MORRIS_CONVERGENCE_BOOTSTRAPS = 200
MORRIS_CONFIDENCE_LEVEL = 0.95
MORRIS_NOISE_MULTIPLIER = 2.0
MORRIS_BOOTSTRAP_SEED = None  # None: use SEED.

# Independent PCE validation thresholds (existing defaults).
PCE_MIN_R2 = 0.90
PCE_MAX_NOISE_MSE = 2.0
PCE_MIN_SIGNAL_TO_NOISE = 3.0
NATIVE_PLOT_SEED = None  # None: use SEED.
HIGH_WORKERS = 1  # Spawned processes for PCE transport; WORKERS is for Morris.
