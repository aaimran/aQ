# Plan: Import `response="anelastic"` (attenuation/Q) from `waveqlab3d_Q_a` into `waveqlab3d_main_Q`

## Goal
Bring the attenuation (anelastic/Q) capability implemented in `waveqlab3d_Q_a` into `waveqlab3d_main_Q`, so that setting `response = "anelastic"` works for:
- `fd_type = "traditional"`
- `fd_type = "upwind"`
- `fd_type = "upwind_drp"`

Scope: focus on the Fortran implementation under `src/` (plus minimal input/README updates needed to run it).

---

## What differs today (src-level delta)
A `src/` diff between the two folders shows these files differ:
- `src/datatypes.f90`
- `src/material.f90`
- `src/block.f90`
- `src/fields.f90`
- `src/RHS_Interior.f90`
- `src/elastic.f90`
- `src/domain.f90`
- `src/moment_tensor.f90` (format/indent only)
- `src/seismogram.f90` (station-list and naming changes; appears unrelated to attenuation)

### Attenuation-relevant deltas in `waveqlab3d_Q_a`
1) **New anelastic state stored in `block_material`** (`src/datatypes.f90`)
- `M%anelastic` flag
- GSLS parameter scalars: `c`, `weight_exp`, `fref` (and `dt` placeholder)
- Per-point fields: `Qp_inv`, `Qs_inv`
- Mechanism arrays (N=4): `tau(4)`, `weight(4)`
- Memory variables for 6 stresses: `eta4..eta9` and their rates `Deta4..Deta9`

2) **Material initialization for anelastic** (`src/material.f90`)
- New subroutine `init_anelastic_properties(M, G, infile)` which:
  - reads namelist `&anelastic_list` (`c`, `weight_exp`, `fref`)
  - allocates/initializes `Qp_inv`, `Qs_inv`, `eta*`, `Deta*`
  - computes `Qp_inv`/`Qs_inv` and sets `tau/weight`
  - adjusts moduli to “unrelaxed” values (updates `M%M(:,:,:,1:2)`)

3) **Block init hook** (`src/block.f90`)
- When `trim(response) == 'anelastic'`, calls `init_anelastic_properties(...)`.

4) **RHS contribution and memory-variable ODEs** (`src/RHS_Interior.f90`)
- In `RHS_near_boundaries`, when `M%anelastic`:
  - subtracts `sum(eta*)` from stress rates
  - updates `M%Deta*` using GSLS evolution equations
- Requires `M` to be `intent(inout)` and propagated from `elastic`.

5) **Time integration of memory variables** (`src/fields.f90`)
- `scale_fields(...)` also scales `M%Deta*`
- `update_fields(...)` advances `M%eta*` using `dt*M%Deta*`

6) **Response parsing/validation** (`src/domain.f90`)
- `waveqlab3d_main_Q` currently validates `response` and only allows `elastic, plastic`.
- `waveqlab3d_Q_a` removed the validation and accepts `anelastic` implicitly.

### Non-attenuation deltas we should NOT blindly port
- `src/seismogram.f90` and the removal of `station_xyz_index` in `src/datatypes.f90` look like an unrelated behavior change (station list parsing + filename scheme). For the attenuation import, we should keep `waveqlab3d_main_Q` behavior unless you explicitly want the seismogram changes.

---

## Important correctness check (potential gap in Q_a)
In `waveqlab3d_Q_a`, the anelastic terms are implemented inside `RHS_near_boundaries`, which *skips* interior points via:
- `if ( ... interior ... ) cycle`

At the same time, the optimized interior RHS (`RHS_Center` → `JJU_x*_interior_*` routines in `src/JU_xJU_yJU_z6.f90`) does **not** appear to include any anelastic/memory-variable logic (no `eta*`, `Qp_inv`, `Qs_inv`, or `anelastic` references).

This implies a risk that, as written, attenuation is only applied near boundaries/PML zones and not throughout the interior.

### Recommended resolution options
- **Option A (best performance):** extend the optimized interior kernels (`JJU_x*_interior_*`) to include the same anelastic stress subtractions and `Deta*` updates.
- **Option B (simplest/safer):** when `response == 'anelastic'`, skip `RHS_Center` and compute the full-domain RHS using the generic looping path (i.e., a modified `RHS_near_boundaries` that does not `cycle` interior points when anelastic). This is slower but likely correct and easier to validate.

The import work should include one of these so `anelastic` is applied everywhere.

---

## Implementation plan (port into `waveqlab3d_main_Q`)

### Step 1 — Add anelastic fields to `block_material` (minimal, keep seismogram fields)
File: `src/datatypes.f90`
- Add the anelastic members from `Q_a` into `type(block_material)`:
  - `anelastic`, `c`, `weight_exp`, `fref`, `Qp_inv`, `Qs_inv`, `tau(4)`, `weight(4)`, `eta4..eta9`, `Deta4..Deta9`
- **Do not remove** existing `seismogram_type` fields like `station_xyz_index` unless you want that feature change.

### Step 2 — Port `init_anelastic_properties` into the material module
File: `src/material.f90`
- Copy `init_anelastic_properties(M, G, infile)` from `Q_a`.
- Ensure it uses the same allocators (`allocate_array_body`) and matches `block_material` member names.
- Decide how to handle missing `&anelastic_list` in input:
  - Either require it (current `Q_a` behavior stops on read error)
  - Or provide a graceful default if the namelist isn’t present

### Step 3 — Hook anelastic init in block setup
File: `src/block.f90`
- Update the `use material, only: ...` line to include `init_anelastic_properties`.
- In `init_block(...)`, add:
  - `if (trim(response) == 'anelastic') call init_anelastic_properties(B%M, B%G, infile)`
  - else ensure `B%M%anelastic = .false.`

### Step 4 — Allow `response="anelastic"` in domain parsing
File: `src/domain.f90`
- Keep `waveqlab3d_main_Q`’s MPI-safe validation, but expand supported responses to:
  - `elastic`, `plastic`, `anelastic`
- Keep the normalize/broadcast logic; do **not** regress to accepting arbitrary strings.

### Step 5 — Integrate memory-variable evolution into time stepping
Files:
- `src/fields.f90`
- `src/elastic.f90`
- `src/RHS_Interior.f90`

Actions:
- Port the `fields.f90` changes to scale and update `M%eta*` / `M%Deta*` when `M%anelastic`.
- Port `RHS_near_boundaries` changes (anelastic stress correction + `Deta*` accumulation) and update intents:
  - `RHS_near_boundaries(..., M)` must take `M` as `intent(inout)`.
  - `set_rates_elastic(..., M, ...)` must pass `M` as `intent(inout)`.

### Step 6 — Fix the “interior attenuation” gap (choose A or B)
- **Option A:** implement the attenuation terms inside each optimized interior kernel used by:
  - traditional: `JJU_x6_interior`
  - upwind: `JJU_x{2..9}_interior_upwind`
  - upwind_drp: `JJU_x{3..7,66,679}_interior_upwind_drp`

  This requires identifying where those kernels compute stress rates and inserting:
  - the `-sum(eta*)` stress corrections
  - the per-mechanism `Deta*` updates

- **Option B:** in `src/elastic.f90`, when `M%anelastic`:
  - skip `RHS_Center(...)` and compute the RHS exclusively via a “full-domain” loop path.
  - simplest is to modify `RHS_near_boundaries` to accept a flag `do_all_points` (or branch on `M%anelastic`) so it does not `cycle` interior points.

Recommendation: start with **Option B** to get a correct baseline quickly; later optimize with Option A if performance is needed.

---

## Validation plan
1) **Build** `waveqlab3d_main_Q` after porting.
2) **Smoke run** a small input for each `fd_type`:
   - `traditional`, `upwind`, `upwind_drp`
   - `response = "anelastic"`
   - include `&anelastic_list c=..., weight_exp=..., fref=... /`
3) Confirm:
   - run does not crash
   - attenuation arrays are allocated (`M%anelastic == .true.`)
   - stresses decay/attenuate relative to `response="elastic"` for the same setup

---

## Deliverables
- Code changes in `waveqlab3d_main_Q/src/` implementing `response="anelastic"` end-to-end.
- (Optional) small input examples in `waveqlab3d_main_Q/inputfile/` showing the new `&anelastic_list` namelist.

---

## Open question (only if you want me to decide for you)
For the initial import, should I implement **Option B** (correctness-first, slower) or go directly for **Option A** (faster, more invasive)?
