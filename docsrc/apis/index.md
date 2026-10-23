# API reference

This section documents every public module, class and function of PyFlange.
It is generated automatically from the docstrings in the source code, so it
always matches the installed version.

| Module | Use it for |
|--------|------------|
| [`pyflange.bolts`](bolts.md) | Bolts, nuts and washers: standard ISO sizes and materials, stiffness and tensile capacity. |
| [`pyflange.flangesegments`](flangesegments.md) | The flange segment models (L-flange and T-flange, AMD1 and ED2): bolt force and bolt moment as functions of the shell force, ULS failure mode check, shell stiffness interpolation. |
| [`pyflange.fatigue`](fatigue.md) | Markov matrices, S-N curves and the bolt fatigue analysis (damage and fatigue life). |
| [`pyflange.stats`](stats.md) | Gap height distributions and random samplers for Monte Carlo simulations. |
| [`pyflange.io`](io.md) | Reading and writing fatigue cases in the [JSON format](../json-format.md), running Monte Carlo simulations of a case. |
| [`pyflange.cli`](cli.md) | The `pyflange` command-line tool. See the [Command line](../cli.md) page. |
| [`pyflange.utils`](utils.md) | Internal helpers: logging, data files, angle reduction. |
| [`pyflange.gap`](gap.md) | Deprecated. Use `pyflange.stats`. |

## How to read this reference

- **Units.** PyFlange does not convert units. Pass all values in one
  consistent system, normally SI (m, N, Pa, rad). See
  [Units and conventions](../guide/units.md).
- **References.** Docstrings cite the documents listed on the
  [References](../references.md) page as `ref. [n]`. The numbers are the
  same in every module. For example, `ref. [2], eq. G.36` is equation G.36
  of IEC 61400-6:2020+AMD1:2025.
- **Symbols.** The parameter names follow the symbols of IEC 61400-6
  Annex G where possible. The [Nomenclature](../guide/nomenclature.md) page
  maps the PyFlange names to the IEC symbols.
- **Private members.** Names starting with an underscore are internal and
  are not listed here. They can change without notice.
