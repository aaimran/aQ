# Anelastic Implementation in waveqlab3d_main_seis

## Model: Generalized Standard Linear Solid (GSLS)

The code implements **viscoelastic attenuation via N=4 Zener (SLS) mechanisms** — the standard approach for nearly-constant-Q wave propagation across a frequency band.

---

## 1. Data Structures (`datatypes.f90`)

Each `block_material` carries:

| Field | Shape | Description |
|-------|-------|-------------|
| `Qp_inv`, `Qs_inv` | `(nx,ny,nz)` | $Q_P^{-1}$, $Q_S^{-1}$ at each grid point |
| `tau(4)` | scalar array | Relaxation times $\tau_l$ for 4 mechanisms |
| `weight(4)` | scalar array | Quadrature weights $w_l$ |
| `eta4`–`eta9` | `(nx,ny,nz,4)` | Memory variables for 6 stress components × 4 mechanisms |
| `Deta4`–`Deta9` | same | Rates $\dot{\eta}$ (RHS for time integration) |

The 6 memory variables correspond to stress components: $\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{xy}, \sigma_{xz}, \sigma_{yz}$.

---

## 2. Initialization (`material.f90`)

### Q from velocity ratio

$$Q_S^{-1} = \frac{1}{c\sqrt{\lambda/\mu}}, \qquad Q_P^{-1} = \frac{1}{2}Q_S^{-1}$$

where `c` is a scale factor from the input namelist `&anelastic_list`.

### Relaxation times

Log-spaced over the target frequency band $[0.08, 15]$ Hz:

$$\tau_l = \exp\!\left[\ln\tau_{\min} + \frac{2l-1}{16}(\ln\tau_{\max}-\ln\tau_{\min})\right], \quad l=1,\ldots,4$$

### Pre-tabulated weight sets (selected by `weight_exp`)

| `weight_exp` | $w_1$ | $w_2$ | $w_3$ | $w_4$ |
|---|---|---|---|---|
| `≈ 0.0` | 1.6126 | 0.6255 | 0.6382 | 1.5969 |
| `≈ 0.6` | 0.0336 | 0.6873 | 0.8767 | 1.5202 |

### Unrelaxed moduli correction

Moduli are shifted so the *relaxed* (low-frequency) modulus matches the input velocity:

$$\mu_\text{unrelax} = \frac{\rho v_S^2}{1 - \displaystyle\sum_l \frac{w_l}{(1 + \omega_\text{ref}^2\tau_l^2)\,Q_S^{-1}}}$$

This corrects for the velocity dispersion introduced by the attenuation mechanism.

---

## 3. RHS: Stress Rates (`RHS_Interior.f90`, `JU_xJU_yJU_z6.f90`)

At every interior point, after computing the elastic stress rates, the memory-variable contribution is **subtracted** from $\dot{\sigma}$:

$$\dot{\sigma}_{ij}^\text{total} = \dot{\sigma}_{ij}^\text{elastic} - \sum_{l=1}^{4} \eta_{ij,l}$$

Simultaneously, the memory variable rates are computed.

**Normal stress components** ($\eta_4, \eta_5, \eta_6$ for $\sigma_{xx}, \sigma_{yy}, \sigma_{zz}$):

$$\dot{\eta}_{kk,l} = \frac{1}{\tau_l}\left[ w_l\!\left(2\mu Q_S^{-1}\,\varepsilon_{kk} + \left((\lambda+2\mu)Q_P^{-1} - 2\mu Q_S^{-1}\right)\mathrm{tr}(\varepsilon)\right) - \eta_{kk,l} \right]$$

**Shear components** ($\eta_7, \eta_8, \eta_9$ for $\sigma_{xy}, \sigma_{xz}, \sigma_{yz}$):

$$\dot{\eta}_{ij,l} = \frac{1}{\tau_l}\left[ w_l\,\mu Q_S^{-1}\,\gamma_{ij} - \eta_{ij,l} \right]$$

where $\gamma_{ij} = \partial_i u_j + \partial_j u_i$ are the engineering shear strains from the velocity gradient components `Ux`, `Uy`, `Uz`.

---

## 4. PML Variant (`apply_anelastic_point_pml`)

Inside PML regions the strain rates are modified before computing memory variable rates. PML auxiliary fields `Qx, Qy, Qz` are subtracted:

$$\varepsilon_{xx}^\text{eff} = U_{x,x} - d_x Q_x^{(4)}, \quad \text{etc.}$$

ensuring attenuation couples correctly with the PML absorbing boundary.

A dispatch routine `apply_anelastic_point_dispatch` selects the PML or standard path per grid point.

---

## 5. Time Integration (`fields.f90`)

Memory variables are integrated as part of the **same RK4 loop** as the wavefield:

| RK4 stage action | Code |
|-----------------|------|
| Scale accumulators | `Deta* = A * Deta*` |
| Update state | `eta* = eta* + dt * Deta*` |

The `eta*` arrays carry state across steps; `Deta*` are accumulated during each RK4 stage and reset.

---

## 6. Input Namelist

Activated when `response = 'anelastic'` in the block input file.

```fortran
&anelastic_list
  c          = 1.0      ! Q scaling: Qs = c * Vp/Vs
  weight_exp = 0.0      ! selects weight set (0.0 or 0.6)
  fref       = 1.0      ! reference frequency (Hz) for unrelaxed modulus correction
/
```

---

## 7. Key Design Choice: No MPI Ghost Exchange for Memory Variables

Memory variables are **point-local** — their evolution equation contains no spatial derivatives of $\eta$. Only the stress/velocity fields require halo (ghost) exchanges. This significantly simplifies the parallel implementation: no additional MPI communication is needed when attenuation is enabled.
