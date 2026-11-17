# Engineering references, extracted parameters and explicit assumptions

Research date: 2026-11-06. This is a DEMO-inspired neutronics model, not an
approved DEMO design or a reconstruction of released engineering CAD. References
from different design years must not be interpreted as one consistent baseline.
Crossref bibliographic metadata is saved in `results/reference_metadata.json`.

## Primary references

1. **R. Kembleton, J. Morris, M. Siccinio, F. Maviglia (2022), “EU-DEMO design
   space exploration and design drivers,” Fusion Engineering and Design.**
   [DOI](https://doi.org/10.1016/j.fusengdes.2022.113080).
   Extracted: example baseline R=9.0 m, a=2.9 m, elongation 1.65,
   triangularity 0.33, fusion power 2000 MW. Radial-build discussion gives
   plasma/first-wall gap 225 mm, IB/OB blanket 755/982 mm, VV total 600 mm,
   thermal shield 50 mm. These build values motivate default parameters.
   Confidence: high for the cited baseline, moderate when combined with [2].

2. **V. Corato, C. Vorpahl, K. Sedlak, V. A. Anvar et al. (2022), “The DEMO
   magnet system – Status and future challenges,” Fusion Engineering and Design.**
   [DOI](https://doi.org/10.1016/j.fusengdes.2021.112971).
   Extracted: R=9.1 m, a=2.9 m, 16 TF coils, approximately 19 m TF overall
   height, 2 GW fusion power; low-temperature Nb3Sn concepts. Model adopts R=9.1 m
   with shaping from [1]. PF/CS positions and detailed pack dimensions below
   are assumptions, not digitized coordinates from this paper.
   Confidence: high for machine scale/count, low for assumed detailed geometry.

3. **F. Hernández, P. Pereslavtsev, Q. Kang, P. Norajitra, B. Kiss, G. Nádasi,
   O. Bitz (2017), “A new HCPB breeding blanket for the EU DEMO: Evolution,
   rationale and preliminary performances,” Fusion Engineering and Design 124,
   882–886.** [DOI](https://doi.org/10.1016/j.fusengdes.2017.02.008),
   [KIT record/full text](https://publikationen.bibliothek.kit.edu/1000067338).
   Extracted: design evolution and modular breeder-unit/cooling architecture.
   This project deliberately uses an earlier sandwich-like Li-bearing ceramic/
   Be bed architecture rather than claiming to implement current pin designs.
   Confidence: high for concept, moderate for the adapted architecture.

4. **I. Palermo, F. A. Hernández, P. Pereslavtsev, D. Rapisarda, G. Zhou (2022),
   “Shielding Design Optimization of the Helium-Cooled Pebble Bed Breeding
   Blanket for the EU DEMO Fusion Reactor,” Energies 15, 5734.**
   [DOI/full text](https://doi.org/10.3390/en15155734).
   Extracted: representative 20 mm FW, 30 mm backplate; current breeder
   mixture includes Li4SiO4 and Li2TiO3, with 60% Li6 enrichment. These support
   FW steel/coolant envelope 2 cm, backplate 3 cm and enrichment 0.6; they do not establish
   Li2O engineering qualification or the assumed internal bed fractions.
   Confidence: high for the reported values; low for Li2O substitution performance.

5. **C. Bachmann, L. Ciupinski, C. Gliss, T. Franke et al. (2022), “Containment
   structures and port configurations,” Fusion Engineering and Design 174, 112966.**
   [DOI](https://doi.org/10.1016/j.fusengdes.2021.112966),
   [open manuscript](https://iris.enea.it/retrieve/408f0e0f-43ad-403c-9ad7-577e93e4eeda/Containment%20structures%20and%20port%20configurations%20_%20Elsevier%20Enhanced%20Reader.pdf).
   Extracted: double-shell austenitic vessel and water/steel shielding;
   inner-shell revisions include 32/36 mm IB/OB. Model uses a uniform 36 mm
   inner shell, assumed 60 mm outer shell and filler to maintain the 600 mm
   total build. Actual port/support geometry is not reconstructed.
   Confidence: moderate for this simplified radial build.

6. **EUROfusion WPHCD-CP(16)15081 submitted manuscript (2016).**
   [Manuscript](https://scipub.euro-fusion.org/wp-content/uploads/eurofusion/WPHCDCP16_15081_submitted.pdf).
   Extracted: HCPB sandwich beds, steel cooling plates, BSS manifolds,
   helium coolant 8 MPa and 300–500°C; purge gas distinct from coolant.
   Full author metadata was not independently recovered from the accessible
   snippet; do not use this identifier as a complete bibliography citation.
   Model ideal-gas helium density uses 8 MPa, 673.15 K; purge pressure 0.2 MPa
   is an assumption. Confidence: moderate for adapted coolant conditions.

7. **G. Stornelli, A. Di Schino, R. Montanari, C. Testani, A. Varone,
   S. Mancini (2021), “Work-Hardening Behavior of Cold Rolled EUROFER97 Steel
   for Nuclear Fusion Applications,” Materials Proceedings 3, 21.**
   [DOI](https://doi.org/10.3390/IEC2M-09242),
   [composition table](https://www.mdpi.com/2673-4605/3/1/21).
   Extracted representative heat (wt%): Cr 8.93, C 0.12, Mn 0.47, V 0.20,
   W 1.07, Ta 0.14, Ti 0.009, N 0.018, small reported trace elements, Fe balance.
   Model omits unspecified upper-bound impurities rather than inventing them.
   Density 7.8 g/cm3 is an engineering assumption. Confidence: moderate.

8. **Particle Data Group, Atomic and Nuclear Properties: lithium oxide.**
   [PDG table](https://pdgprod.lbl.gov/~gtest/2026/AtomicNuclearProperties/HTML/lithium_oxide.html).
   Extracted solid density 2.013 g/cm3 and Li2O stoichiometry. This is a reference
   density, not a validated high-temperature swelling/thermal-expansion law.
   Confidence: high for reference density, low for irradiated hot-bed properties.

9. **J. H. You et al., “Structural lifetime assessment for the DEMO divertor
   targets: Design-by-analysis” (open manuscript).**
   [Manuscript](https://pure.mpg.de/pubman/item/item_3277677_4/component/file_3278542/You_Structural.pdf).
   Extracted: W armour, Cu/CuCrZr heat-sink/cooling concept; target armour
   alternatives 5–8 mm. Model uses 8 mm W with simplified inclined target
   strips; target/cassette positions and cooling-strip geometry are assumptions.
   The full journal metadata should be checked before formal publication.
   Confidence: moderate for materials/armour, low for placement.

10. **OpenMC developers, installed source and official CSG API documentation.**
    [Polygon](https://docs.openmc.org/en/stable/pythonapi/generated/openmc.model.Polygon.html),
    [CSG methods](https://docs.openmc.org/en/stable/methods/geometry.html).
    Audited local `openmc/model/surface_composite.py`: `Polygon(basis='rz')`
    decomposes polygons into native planes, cylinders and one-sided cones.
    `openmc/source.py` and `tests/unit_tests/test_source_tokamak.py` define the
    installed Miller-shaped source. No arbitrary elliptical-torus primitive
    or DAGMC transport is assumed. Confidence: high for installed capabilities;
    runtime tests, not documentation alone, establish usable behaviour.

## Assumption register

| Parameter/design choice | Default | Status/confidence |
|---|---:|---|
| Major/minor radius | 910/290 cm | Magnet baseline [2]; high |
| Elongation/triangularity | 1.65/0.33 | [1], transplanted to R=9.1 m; moderate |
| Plasma contour | 32 Miller vertices | Numerical faceting assumption; refinement required |
| FW gap | 22.5 cm | [1]; high for that baseline |
| IB/OB blanket | 75.5/98.2 cm | [1]; high for that baseline |
| FW armour | 0.2 cm W | Published DEMO concept value; moderate across variants |
| FW steel/cooling split | 0.8/1.0/0.2 cm | Assumed split within 2 cm steel/coolant envelope |
| Blanket segmentation | 16 TF sectors × 3 toroidal × 16 poloidal bins | Assumption; lower bins truncated by divertor opening |
| Toroidal/poloidal gaps | 1 cm / 1 cm | Assumptions; vertical IB/OB stacks, poloidal sidewall 0.5 cm |
| Active radial cycles | 3 | Assumed design study parameter |
| Active material fractions | breeder 0.48, multiplier 0.45, cooling 0.07 | Assumed fractions, sum one |
| Bed packing | 0.64 | Explicit typical engineering assumption; not material solid density |
| Rib fraction | 0.03 | Assumed equivalent structural volume |
| Manifold/backplate | 12/3 cm | Manifold assumed; backplate motivated by [4] |
| Breeder | Li2O, 60% Li6 | Deliberate departure from HCPB reference ceramic |
| Multiplier | Be, alternative Be12Ti | Be sandwich baseline; alternative density assumed 2.28 g/cm3 |
| TF winding pack/case | 75 cm build, 70 cm toroidal width | Assumption adjusted to DEMO scale; not stress-optimized |
| TF D shape | minor 590 cm, elongation 1.45, triangularity 0.15 | Geometric engineering envelope assumption |
| PF coils | Six individual rings at configured R/Z | Assumed positions outside TF envelopes |
| CS | Six discrete modules, R=65–165 cm, total height 1080 cm | Assumption; not a magnetic design |
| Winding material | SS/Cu/Nb3Sn/glass/He = 35/35/10/10/10 vol% | Explicit homogenization assumption |
| Insulation | Silica-glass equivalent | Reduced chemistry; no validated resin radiation model |
| VV | 3.6 cm inner, 50.4 cm filler, 6 cm outer | Total from [1], partial split adapted [5] |
| VV filler | 60 vol% SS, 40% water | Assumed shielding/cooling mix |
| Divertor opening | Blanket omitted below Z=−400 cm | Geometric single-null approximation; not an equilibrium |
| Divertor target positions | Coordinates in `geometry/divertor.py` | Assumed V-like inclined surfaces; verify after dimensional changes |
| Water | 0.7 g/cm3, no thermal S(a,b) | Approximate hot-water neutronics; not thermal-hydraulic calculation |
| Magnet temperature | transport 293.6 K, operating label 4.5 K | Explicit XS-temperature approximation; not a cryogenic calculation |
| Nuclear-data temperature | nearest, configured tolerance | Record selected native temperatures; interpolation/convergence remains future work |
| Fusion source | 14.1 MeV, isotropic, uniform plasma volume | Standard DT-energy approximation; no beam angular anisotropy |

## Li2O substitution: interpretation

Li2O replaces the reference lithium silicate/titanate ceramic; the resulting model
must not be labelled the evaluated EUROfusion HCPB reference. At fixed bed volume,
solid density, stoichiometry and enrichment determine lithium atom inventory.
The bed contains explicit solid-volume fraction and interstitial purge helium;
its effective density is calculated by volume mixing, not by using solid density
for the whole bed. Oxygen natural isotopes are retained where the library supports
them. Enrichment changes isotopic atom fractions, not an arbitrary cross-section
factor. Reference solid density is kept configurable; thermal expansion is not
silently invented. Tritium release, hygroscopic behaviour, chemical compatibility,
thermal conductivity, pebble mechanics, swelling and thermo-mechanical safety
require separate models and qualification. A neutronic TBR is not proof of
breeder thermal feasibility or fuel-cycle self-sufficiency.
