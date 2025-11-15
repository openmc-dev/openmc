import openmc

sp = openmc.StatePoint('statepoint.150.h5')

print("\nResults:")
print("=" * 70)
print(f"k-effective:              {sp.keff}")
print(f"k-prompt:                 {sp.k_prompt}")
print(f"Beta-effective:           {sp.beta_eff}")
print(f"Prompt generation time:   {sp.prompt_gen_time} seconds")
print(f"Alpha (k-based):          {sp.alpha_k_based} 1/s")
print(f"Alpha (rate-based):       {sp.alpha_rate_based} 1/s")

# Expected results for Godiva (from Cullen et al. 2003):
print("\nExpected values for Godiva (UCRL-TR-201506):")
print("  k_eff ≈ 1.0 (near critical)")
print("  beta_eff ≈ 0.0065-0.0070 (0.65-0.70%)")
print("  k_prompt ≈ 0.993-0.994")
print("  This problem demonstrates static vs dynamic criticality concepts.")
