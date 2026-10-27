# Nomenclature

PyFlange names follow the symbols of IEC 61400-6 Annex G (ref. [2]) where
possible. This page maps the parameters and properties of the flange segment
classes to the symbols of the standard.

## Geometry

=== "L-flange"

    <figure>
    --8<-- "assets/img/lflange-geometry.svg"
    <figcaption>Cross-section of an L-flange segment (one bolt). The flange is on the inside of the tower.</figcaption>
    </figure>

=== "T-flange"

    <figure>
    --8<-- "assets/img/tflange-geometry.svg"
    <figcaption>Cross-section of a T-flange segment (one bolt pair). For T-flanges, a is measured from the bolt axis to the edge of the flange.</figcaption>
    </figure>

The drawings are schematic: proportions are not to scale, and the fillet
radius `r`, washers and gap are not shown. Z is the shell force acting on the
segment, positive in tension.

## Flange segment parameters

| PyFlange | IEC 61400-6 | Unit | Description |
|----------|-------------|------|-------------|
| `a` | a | m | L-flange: distance from the inner face of the flange to the bolt axis. T-flange: distance from the flange edge to the bolt axis. |
| `b` | b | m | Distance from the bolt axis to the center line of the shell. |
| `s` | s | m | Shell (tower wall) thickness. |
| `s_ratio` | – | – | Thickness of the lower shell divided by `s`. The mean thickness is used where the standard requires it. |
| `t` | t | m | Flange thickness. |
| `R` | R, D/2 | m | Outer radius of the shell. |
| `central_angle` | – | rad | Angle of the segment: 2π / n<sub>bolts</sub> (L-flange) or 2π / n<sub>pairs</sub> (T-flange). The segment width at the shell mid-surface is c = central_angle · (R − s/2). |
| `r` | r | m | Fillet radius between flange and shell. |
| `Zg` | Z<sub>1</sub>, Z<sub>dw</sub> | N | Shell force of the segment at rest (dead weight); negative for compression. |
| `Fv` | F<sub>p,C</sub> (before losses) | N | Bolt preload before settlement and plastic losses. |
| `PPL` | ΔF<sub>pl</sub> / F<sub>p</sub> | – | Plastic preload loss, as a fraction of `Fv` (ref. [2], eq. (24)). |
| `DFT` | DFT<sub>sbw</sub> | m | Coating thickness under the washers, for the settlement loss (ref. [2], eq. (28)). |
| `settlement_factor` | – | – | 1.0, or 0.5 for tension tightening with re-tightening (ref. [2], 6.7.5.1.2). |
| `Do` | d<sub>B</sub> | m | Bolt hole diameter. |
| `bolt`, `nut`, `washer` | – | – | Bolt, nut and washer objects (`pyflange.bolts`). The washer outer diameter, or the nut bearing diameter without washer, is D<sub>w</sub>. |
| `gap.height` | k, k<sub>design</sub> | m | Gap height. |
| `gap.angle` | α<sub>gap</sub> | rad | Gap angle. The gap length is L<sub>gap</sub> = R · α<sub>gap</sub>. |
| `gap.shape_factor` | – | – | Scales the bolt force curve for gaps that are not sinusoidal (1.0 = sinusoidal). |
| `tilt_angle` | α<sub>incl</sub> | rad | Inclination of the flange contact surface (ref. [2], G.2.4). |
| `E`, `G` | E, G | Pa | Elastic and shear modulus of the flange. Poisson's ratio is derived as ν = E/(2G) − 1. |
| `k_shell` | k<sub>shell,ini</sub> | N/m/m | Initial shell stiffness; if `None`, eq. G.43 is used. |

## Results and intermediate quantities

| PyFlange | IEC 61400-6 | Description |
|----------|-------------|-------------|
| `bolt_axial_force(Z)` | F<sub>S</sub>(Z) | Bolt axial force. |
| `bolt_bending_moment(Z)` | M<sub>S</sub>(Z) | Bolt bending moment. |
| `bolt_preload` | F'<sub>p,C</sub> | Preload after losses (ref. [2], eq. (23) without the COV term). |
| `shell_force_at_rest`, `bolt_force_at_rest` | Z<sub>1</sub>, F<sub>S,1</sub> | Point 1. |
| `shell_force_at_tensile_ULS`, `bolt_force_at_tensile_ULS` | Z<sub>2</sub>, F<sub>S,2</sub> | Point 2. |
| `shell_force_at_small_displacement`, `bolt_force_at_small_displacement` | Z<sub>3</sub>, F<sub>S,3</sub> | Point 3. |
| `shell_force_at_closed_gap` | Z<sub>close</sub> | Point 4: shell force that closes the gap. |
| `bolt_moment_at_rest`, `..._at_tensile_ULS`, `..._at_small_displacement` | M<sub>S,1</sub>, M<sub>S,2</sub>, M<sub>S,3</sub> | Bolt moments at the points. |
| `failure_mode(...)` | F<sub>u,A</sub>, F<sub>u,B</sub>, F<sub>u,D</sub>, F<sub>u,E</sub> | Ultimate shell forces of the failure modes (ref. [2], G.1.2). |

The JSON result files (`pyflange determ`) also report internal quantities,
such as `a_prime` (a'), `Z0` (Z<sub>0</sub>), `DZ_gap` (ΔZ<sub>gap</sub>),
`DZ_gap_inclination`, `DZ_gap_tot`, `stiffness_correction_factor`
(α<sub>k</sub>), `polynomial_initial_slope` (χ<sub>ini,mod</sub>), `gap_k_fac`
(k<sub>fac</sub>), `gap_k_shell` (k<sub>shell,ini</sub>), `gap_k_fl`
(k<sub>fl</sub>), `gap_stiffness_correction_factor` (f<sub>tot</sub>), `a_star`
(a\*) and `I_tg` (I<sub>tg</sub>).

## Fatigue

| PyFlange | IEC 61400-6 | Description |
|----------|-------------|-------------|
| `MarkovMatrix.range`, `.mean`, `.cycles` | ΔM, M<sub>mean</sub>, n | Load bins. |
| `BoltFatigueCurve(diameter, DS_ref, gamma_M)` | d, Δσ<sub>C</sub>, γ<sub>M</sub> | Bolt S-N curve (ref. [2], 6.7.5.4). |
| `BoltFatigueAnalysis.damage` | D | Miner's sum. |
| `BoltFatigueAnalysis.allowable_damage` | D<sub>allow</sub> | Allowable damage (1.0 in ref. [2]). |
| `SMF` | SMF | Stress multiplication factor (model uncertainty, ref. [5], 6.5). |
