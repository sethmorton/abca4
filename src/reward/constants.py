"""Constants for reward optimization slice."""

# Greedy Selection
GREEDY_COVERAGE_BONUS = 0.05

# Cross-Entropy Method (CEM)
CEM_ITERATIONS = 300
CEM_POPULATION_SIZE = 120
CEM_ELITE_FRACTION = 0.2
CEM_LAPLACE_SMOOTHING = 1e-6

# Optimization Targets
DEFAULT_QUANTILE_TARGET = 0.7
MAX_TARGET_SCORE = 0.9
