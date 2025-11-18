# Alpha Eigenvalue Calculation in OpenMC

## Overview

The alpha eigenvalue (α) represents the time rate of change of the neutron population in a nuclear system. OpenMC calculates alpha using kinetics parameters derived from prompt neutron chain tallies during normal eigenvalue iterations.

**Formula:**
```
α = (k_prompt - 1) / l_prompt
```

where:
- **k_prompt** = prompt neutron multiplication factor (excluding delayed neutron chains)
- **l_prompt** = prompt neutron lifetime (average time from birth to ANY absorption)

---

## Physical Meaning

The alpha eigenvalue describes population dynamics:
- **α > 0**: Supercritical (k_prompt > 1) → population growing exponentially
- **α = 0**: Critical (k_prompt = 1) → population stable
- **α < 0**: Subcritical (k_prompt < 1) → population decaying exponentially

For example, if α = -1.24 × 10⁶ s⁻¹, the prompt neutron population decays at a rate of 1.24 million per second.

---

## Implementation

### 1. Tally Setup

During initialization, if `calculate_alpha = True`, OpenMC creates an internal tally with prompt-chain scores:

```cpp
void setup_kinetics_tallies()
{
  auto* tally = Tally::create();
  tally->set_writable(false);  // Internal use only

  vector<std::string> scores;
  scores.push_back("prompt-chain-gen-time-num");    // Numerator: Σ(lifetime × weight)
  scores.push_back("prompt-chain-gen-time-denom");  // Denominator: Σ(weight)
  scores.push_back("prompt-chain-nu-fission-rate"); // For diagnostics
  scores.push_back("prompt-chain-absorption-rate"); // For diagnostics
  scores.push_back("prompt-chain-leakage-rate");    // For diagnostics
  scores.push_back("prompt-chain-population");      // For diagnostics

  tally->set_scores(scores);
  tally->set_filters({});  // Tally over entire geometry
}
```

**Key concept**: Prompt-chain scores only accumulate for neutrons born from prompt fission (excluding delayed neutron chains). This allows separate tracking of prompt vs total behavior.

### 2. Calculation During Active Batches

For each active generation, `calculate_kinetics_parameters()` computes:

#### Step 1: Prompt Neutron Lifetime

```cpp
// Extract tally results (accumulated over active batches)
double gen_time_num = results(0, 0, SUM) / n_active;
double gen_time_denom = results(0, 1, SUM) / n_active;

// Calculate prompt neutron lifetime
l_prompt = gen_time_num / gen_time_denom
```

**Physical meaning**:
- Numerator: Sum of (neutron lifetime × neutron weight) for all prompt-chain neutrons
- Denominator: Sum of neutron weights
- Result: Weighted average lifetime from birth to ANY absorption (capture or fission)

**Important distinction**:
- **Prompt neutron lifetime (l_prompt)**: Time from birth to ANY absorption
- **Prompt generation time (Λ_prompt)**: Time from birth to FISSION only
- We use l_prompt, not Λ_prompt

#### Step 2: Prompt k-effective

k_prompt is calculated separately during the simulation by tracking prompt-chain fission banks. It represents the multiplication factor when delayed neutrons are excluded.

#### Step 3: Beta-Effective

```cpp
// β_eff = (k_eff - k_prompt) / k_eff
beta_eff = (keff - keff_prompt) / keff
```

This represents the fraction of fissions producing delayed neutrons.

#### Step 4: Alpha Eigenvalue

```cpp
// α = (k_prompt - 1) / l_prompt
alpha = (keff_prompt - 1.0) / prompt_gen_time
```

### 3. Uncertainty Propagation

Standard deviations are calculated using error propagation for derived quantities:

**For l_prompt:**
```
σ_l² ≈ (∂l/∂num)² σ_num² + (∂l/∂denom)² σ_denom²
```

**For α:**
```
σ_α² ≈ (1/l)² σ_k² + ((k-1)/l²)² σ_l²
```

---

## Tally Scores Explained

### Active Scores (Used in Calculation)

| Index | Score | Description |
|-------|-------|-------------|
| 0 | `prompt-chain-gen-time-num` | Σ(lifetime × weight) - numerator for l_prompt |
| 1 | `prompt-chain-gen-time-denom` | Σ(weight) - denominator for l_prompt |

### Diagnostic Scores (Not Currently Used)

| Index | Score | Description |
|-------|-------|-------------|
| 2 | `prompt-chain-nu-fission-rate` | Fission production rate |
| 3 | `prompt-chain-absorption-rate` | Total absorption rate |
| 4 | `prompt-chain-leakage-rate` | Leakage rate |
| 5 | `prompt-chain-population` | Prompt neutron population |

These diagnostic scores could be used for alternative alpha calculation methods (rate-based: α = (Production - Removal) / Population) but are not currently utilized.

---

## Example Calculation

For a typical fast system (Godiva):

**Given:**
- k_prompt = 0.993
- l_prompt = 5.66 × 10⁻⁶ seconds (5.66 microseconds)

**Calculate:**
```
α = (0.993 - 1.0) / (5.66 × 10⁻⁶)
α = -0.007 / (5.66 × 10⁻⁶)
α = -1.237 × 10⁶ s⁻¹
```

This means the prompt neutron population decays at 1.237 million per second, requiring delayed neutrons to sustain criticality.

---

## Output

OpenMC prints alpha results in the summary:

```
 Delayed Neutron Kinetics Parameters:
 k-effective (Collision)     = 1.00000 +/- 0.00050
 k-prompt (Collision)        = 0.99300 +/- 0.00045
 Beta-effective              = 0.00700 +/- 0.00010
 Prompt Neutron Lifetime     = 5.66000e-06 +/- 2.50000e-08 seconds
 Alpha Eigenvalue            = -1.23700e+06 +/- 1.80000e+04 1/seconds
```

The alpha value is also written to statepoint files for post-processing.

---

## Implementation Notes

### Why This Method Works

1. **No population control needed**: Alpha is calculated from k_prompt, which is measured during normal eigenvalue iterations
2. **Same statistics as k_eff**: Alpha gets the same number of samples as k_prompt
3. **Stable for all systems**: Works for deeply subcritical, critical, or supercritical systems
4. **Proper uncertainties**: Error propagation from measured quantities

### Code Location

- **Setup**: `openmc/src/eigenvalue.cpp::setup_kinetics_tallies()`
- **Calculation**: `openmc/src/eigenvalue.cpp::calculate_kinetics_parameters()`
- **Output**: `openmc/src/output.cpp::print_results()`

### Variable Naming Note

The variable `prompt_gen_time` in the code stores the **prompt neutron lifetime** (l_prompt), NOT the prompt generation time (Λ_prompt). These are different physical quantities:
- **l_prompt**: Average time from birth to ANY absorption (what we calculate)
- **Λ_prompt**: Average time from birth to FISSION (different concept)

---

## References

The relationship α = (k_prompt - 1) / l_prompt is derived from reactor kinetics theory and represents the fundamental mode decay constant for the prompt neutron population.
