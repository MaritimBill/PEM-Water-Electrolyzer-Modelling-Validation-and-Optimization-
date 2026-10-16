
---

## File Placement

| File | Location | Purpose |
|---|---|---|
| `comparison.m` | `PlatEMO/` (same level as `Algorithms/`, `Problems/`) | Driver script — runs the 19-algorithm comparison and the post-hoc Simulink validation |
| `PEM_Problem_Scaled.m` | `PlatEMO/Problems/Your Problem/` | The 3-objective constrained PEM design problem |
| `NSGAII_AdaptiveSearch.m` | `PlatEMO/Sims/NSGAII_AdaptiveSearch` | Proposed NSGA-II variant (adaptive ops + local search + diversity injection) |
| `EnvironmentalSelection.m` | `PlatEMO/Sims/NSGAII_AdaptiveSearch` | Custom environmental selection used by the variant |
| `build_PEM_Simulink.m` | `PlatEMO/Sims/` | Builds `PEM_Plant.slx` from scratch if missing |
| `PEM_Simulink_Wrapper.m` | `PlatEMO/Sims/` | Accelerator + Fast-Restart wrapper |
| `test_simulink_wrapper.m` | `PlatEMO/Sims/` | Standalone check: analytical model vs. Simulink |
| `PEM_Plant.slx` | `PlatEMO/Sims/` | Simulink digital twin |

**Rules of thumb:**

- Everything stock PlatEMO stays where it is. Do not copy algorithm
  files, `ALGORITHM.m`, `PROBLEM.m`, or any PlatEMO internals into `Sims/`.
- `comparison.m` sits at the PlatEMO root so `fileparts(mfilename('fullpath'))`
  resolves to the PlatEMO root directly.
- `Sims/` is a subfolder of PlatEMO, added to the path by `genpath`. Its
  contents are the only files unique to this project.

---

## Setup

1. **Install PlatEMO v4.16 (2026_07_19).** Extract it wherever you like;
   note the path to the folder that contains `platemo.m`.

2. **Copy the files from this repo** into the PlatEMO tree per the
   placement table above. In particular:
   - `comparison.m` → PlatEMO root
   - `PEM_Problem_Scaled.m` → `PlatEMO/Problems/Your Problem/`
   - everything else → `PlatEMO/Sims/`

3. **Open MATLAB R2023b(or one you have)** and `cd` into the PlatEMO root:

   ```matlab
   cd('<path-to>/PlatEMO')
