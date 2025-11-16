import openmc
import warnings

# Suppress the cross sections path warning
warnings.filterwarnings('ignore', message="Path.*does not exist", category=UserWarning)

sp = openmc.StatePoint('statepoint.150.h5')

# Convert units for fast neutron systems
# StatePoint stores values in SI units (seconds), convert to microseconds for display
prompt_gen_time_us = sp.prompt_gen_time.nominal_value * 1e6 if sp.prompt_gen_time else 0.0
prompt_gen_time_std_us = sp.prompt_gen_time.std_dev * 1e6 if sp.prompt_gen_time else 0.0

alpha_k_based_us = sp.alpha_k_based.nominal_value / 1e6 if sp.alpha_k_based else 0.0
alpha_k_based_std_us = sp.alpha_k_based.std_dev / 1e6 if sp.alpha_k_based else 0.0

print("\nResults:")
print("=" * 70)
print(f"k-effective:              {sp.keff}")
print(f"k-prompt:                 {sp.k_prompt}")
print(f"Beta-effective:           {sp.beta_eff}")
print(f"Prompt generation time:   {prompt_gen_time_us:.6e}+/-{prompt_gen_time_std_us:.6e} us")
print(f"Alpha (k-based):          {alpha_k_based_us:.6e}+/-{alpha_k_based_std_us:.6e} gen/us")

# Expected results for Godiva (from Cullen et al. 2003):
print("\nExpected values for Godiva (UCRL-TR-201506):")
print("  k_eff ≈ 1.0 (near critical)")
print("  beta_eff ≈ 0.0065-0.0070 (0.65-0.70%)")
print("  k_prompt ≈ 0.993-0.994")
print("  Prompt generation time ≈ 5-10 ns (5e-3 to 1e-2 us for fast systems)")
print("  This problem demonstrates static vs dynamic criticality concepts.")
