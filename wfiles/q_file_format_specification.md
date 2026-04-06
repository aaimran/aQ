# File Format Specifications for Withers Q Models

## Overview

This document specifies the binary file formats for reading spatially-varying Q and gamma fields when using `Q0_type='file'` in the Withers implementation.

---

## 1. Q Model Files (Qs0_file, Qp0_file)

### Format: Direct Access Binary (Fortran Unformatted)

**Purpose:** Store spatially-varying quality factor Q₀ values

**File Extension:** `.bin` or `.dat`

**Byte Order:** Native (machine-dependent, typically little-endian)

**Data Type:** `real(kind=wp)` where `wp` is working precision (typically `real*8`, 64-bit double precision)

### Structure

```
┌─────────────────────────────────────────────────┐
│ Header (optional, can be skipped)              │
├─────────────────────────────────────────────────┤
│ 3D Array: Q(i,j,k)                             │
│   Dimensions: [nq, nr, ns]                     │
│   Order: Fortran column-major (k fastest)      │
│   Values: Q₀ at each grid point               │
└─────────────────────────────────────────────────┘
```

### Array Ordering

**Fortran column-major (k varies fastest):**
```fortran
! Storage order in file:
! Q(1,1,1), Q(1,1,2), ..., Q(1,1,ns),
! Q(1,2,1), Q(1,2,2), ..., Q(1,2,ns),
! ...
! Q(1,nr,1), ..., Q(1,nr,ns),
! Q(2,1,1), ..., Q(2,nr,ns),
! ...
! Q(nq,nr,ns)
```

### Value Constraints

- **Minimum Q:** 15.0 (lower bound of Withers method validity)
- **Maximum Q:** 10,000 (essentially elastic)
- **Typical ranges:**
  - Sediments: 10-50
  - Shallow crust: 50-150
  - Deep crust: 200-600
  - Mantle: 600-1000

### Read Example (Fortran)

```fortran
subroutine read_q_file(filename, nq, nr, ns, Q_array)
   use mpi3dio, only: file_distributed, open_file_distributed, &
                      read_file_distributed, close_file_distributed
   use mpi3dcomm, only: cartesian3d_t
   
   character(*), intent(in) :: filename
   integer, intent(in) :: nq, nr, ns
   real(kind=wp), dimension(nq,nr,ns), intent(out) :: Q_array
   
   type(file_distributed) :: fid
   
   ! Open file for parallel reading
   call open_file_distributed(fid, filename, "read", comm, array_desc, pw)
   
   ! Read Q values (assumes no header, starts at byte 1)
   call read_file_distributed(fid, Q_array)
   
   ! Close file
   call close_file_distributed(fid)
   
   ! Validate Q values
   where (Q_array < 15.0_wp) Q_array = 15.0_wp
   where (Q_array > 10000.0_wp) Q_array = 10000.0_wp
   
end subroutine read_q_file
```

### Write Example (Python - for creating test files)

```python
import numpy as np

def write_q_file(filename, Q_array):
    """
    Write Q array to binary file in Fortran-compatible format
    
    Parameters:
    -----------
    filename : str
        Output filename
    Q_array : ndarray
        3D array of Q values, shape (nq, nr, ns)
    """
    # Ensure correct data type (64-bit float, Fortran double precision)
    Q_fortran = np.asfortranarray(Q_array, dtype=np.float64)
    
    # Write as raw binary (no header)
    with open(filename, 'wb') as f:
        f.write(Q_fortran.tobytes(order='F'))
    
    print(f"Wrote Q file: {filename}")
    print(f"  Dimensions: {Q_array.shape}")
    print(f"  Q range: [{Q_array.min():.1f}, {Q_array.max():.1f}]")
    print(f"  File size: {Q_fortran.nbytes / 1024**2:.2f} MB")

# Example: Create simple velocity-based Q model
nq, nr, ns = 200, 150, 100
x = np.linspace(0, 50, nq)  # km
y = np.linspace(0, 40, nr)
z = np.linspace(0, 20, ns)

# Simple depth-dependent Q
Z = z[np.newaxis, np.newaxis, :]
Qs = 50 + 10 * Z  # Q increases with depth

write_q_file('Qs_model.bin', Qs)
```

---

## 2. Gamma File (gamma_file) - Optional

### Format: Direct Access Binary (Fortran Unformatted)

**Purpose:** Store spatially-varying power-law exponent γ

**File Extension:** `.bin` or `.dat`

**Data Type:** `real(kind=wp)` (64-bit double precision)

### Structure

```
┌─────────────────────────────────────────────────┐
│ 3D Array: gamma(i,j,k)                         │
│   Dimensions: [nq, nr, ns]                     │
│   Order: Fortran column-major                   │
│   Values: γ ∈ [0.0, 0.9] at each grid point    │
└─────────────────────────────────────────────────┘
```

### Value Constraints

- **Minimum γ:** 0.0 (frequency-independent Q)
- **Maximum γ:** 0.9 (strong frequency dependence)
- **Typical values:**
  - California: 0.6 - 0.7
  - Eastern US: 0.3 - 0.4
  - High-Q stable regions: 0.1 - 0.3
  - Low-Q tectonic regions: 0.7 - 0.9

### Read/Write Example

**Same as Q files, just different value range**

```python
def write_gamma_file(filename, gamma_array):
    """Write gamma array (same format as Q files)"""
    # Clamp to valid range
    gamma_clamped = np.clip(gamma_array, 0.0, 0.9)
    
    gamma_fortran = np.asfortranarray(gamma_clamped, dtype=np.float64)
    
    with open(filename, 'wb') as f:
        f.write(gamma_fortran.tobytes(order='F'))
    
    print(f"Wrote gamma file: {filename}")
    print(f"  Gamma range: [{gamma_clamped.min():.2f}, {gamma_clamped.max():.2f}]")

# Example: Regionally varying gamma
gamma = np.full((nq, nr, ns), 0.6)  # Base value for California

# Make sedimentary basins have higher gamma
shallow_basin = (z < 2.0)  # Top 2 km
gamma[:, :, shallow_basin] = 0.75

write_gamma_file('gamma_california.bin', gamma)
```

---

## 3. Alternative Format: HDF5 (Future Enhancement)

### Advantages
- Self-describing (metadata included)
- Compression support
- Portable across platforms
- Easy to inspect with tools

### Example Structure

```
file.h5
├── Qs                  # Dataset [nq, nr, ns], float64
│   ├── units: "dimensionless"
│   ├── description: "Shear quality factor Q_s"
│   └── valid_range: [15.0, 10000.0]
├── Qp                  # Dataset [nq, nr, ns], float64
│   ├── units: "dimensionless"
│   └── description: "P-wave quality factor Q_p"
├── gamma               # Dataset [nq, nr, ns], float64 (optional)
│   ├── units: "dimensionless"
│   ├── description: "Power-law exponent"
│   └── valid_range: [0.0, 0.9]
└── metadata
    ├── nq: 200
    ├── nr: 150
    ├── ns: 100
    ├── origin: [0.0, 0.0, 0.0]  # km
    ├── spacing: [0.25, 0.25, 0.2]  # km
    └── coordinate_system: "Cartesian"
```

---

## 4. File Validation and Error Handling

### Validation Checks (Implement in Fortran)

```fortran
subroutine validate_q_file(Q_array, name)
   real(kind=wp), dimension(:,:,:), intent(inout) :: Q_array
   character(*), intent(in) :: name
   
   integer :: nq, nr, ns, n_clamped
   real(kind=wp) :: qmin, qmax
   
   nq = size(Q_array, 1)
   nr = size(Q_array, 2)
   ns = size(Q_array, 3)
   
   qmin = minval(Q_array)
   qmax = maxval(Q_array)
   
   ! Check for invalid values
   n_clamped = 0
   
   where (Q_array < 15.0_wp)
      Q_array = 15.0_wp
      n_clamped = n_clamped + 1
   end where
   
   where (Q_array > 10000.0_wp)
      Q_array = 10000.0_wp
      n_clamped = n_clamped + 1
   end where
   
   ! Report
   if (is_master()) then
      print *, "Loaded Q file: ", trim(name)
      print *, "  Dimensions: ", nq, "x", nr, "x", ns
      print *, "  Q range: [", qmin, ",", qmax, "]"
      if (n_clamped > 0) then
         print *, "  WARNING: Clamped ", n_clamped, " values to [15, 10000]"
      end if
   end if
   
end subroutine validate_q_file
```

---

## 5. Complete Usage Examples

### Example 1: Uniform Q from Velocity (No Files)

```fortran
&withers_list
  gamma = 0.6,
  Q0_type = 'velocity',
  Qs0_velocity_coef = 30.0,  ! Qs = 30 * Vs
  Qp0_over_Qs0 = 2.0,
/
```

### Example 2: Constant Q (No Files)

```fortran
&withers_list
  gamma = 0.0,  ! Constant Q
  Q0_type = 'constant',
  Qs0_constant = 100.0,
  Qp0_over_Qs0 = 2.0,
/
```

### Example 3: Q from Files, Uniform Gamma

```fortran
&withers_list
  gamma = 0.65,  ! Uniform California exponent
  Q0_type = 'file',
  Qs0_file = 'models/scec_cvm_Qs.bin',
  Qp0_file = 'models/scec_cvm_Qp.bin',
  ! gamma_file not specified -> use uniform gamma
/
```

### Example 4: Both Q and Gamma from Files

```fortran
&withers_list
  gamma = 0.6,  ! Fallback value if file read fails
  Q0_type = 'file',
  Qs0_file = 'models/tomography_Qs.bin',
  Qp0_file = 'models/tomography_Qp.bin',
  gamma_file = 'models/gamma_regional.bin',
/
```

---

## 6. Creating Q Models from Velocity Models

### Python Script: Velocity-Based Q

```python
#!/usr/bin/env python3
"""
Create Q model files from velocity model
Using empirical relationships
"""
import numpy as np

def velocity_to_q(vs, formula='withers'):
    """
    Convert Vs to Qs using empirical formulas
    
    Parameters:
    -----------
    vs : ndarray
        S-wave velocity in km/s
    formula : str
        'withers' : Qs = 30 * Vs  (Withers et al., 2015)
        'olsen' : Qs = 13.9 * Vs^2.18  (Olsen et al., 2003)
        'brocher' : Qs = 20 * Vs  (Simple scaling)
    """
    if formula == 'withers':
        return 30.0 * vs
    elif formula == 'olsen':
        return 13.9 * (vs ** 2.18)
    elif formula == 'brocher':
        return 20.0 * vs
    else:
        raise ValueError(f"Unknown formula: {formula}")

def create_q_model(vs_model, vp_model, output_prefix='q_model'):
    """
    Create Qs and Qp binary files from velocity model
    
    Parameters:
    -----------
    vs_model : ndarray
        3D array of Vs in km/s
    vp_model : ndarray
        3D array of Vp in km/s
    output_prefix : str
        Prefix for output files
    """
    # Compute Qs from Vs
    Qs = velocity_to_q(vs_model, formula='withers')
    
    # Compute Qp (typically 2*Qs)
    Qp = 2.0 * Qs
    
    # Clamp to valid range
    Qs = np.clip(Qs, 15.0, 10000.0)
    Qp = np.clip(Qp, 15.0, 10000.0)
    
    # Write files
    write_q_file(f'{output_prefix}_Qs.bin', Qs)
    write_q_file(f'{output_prefix}_Qp.bin', Qp)
    
    return Qs, Qp

# Example usage:
if __name__ == '__main__':
    # Load velocity model (example)
    # vs = np.load('velocity_model_vs.npy')  # Shape: (nq, nr, ns)
    # vp = np.load('velocity_model_vp.npy')
    
    # Or create simple test model
    nq, nr, ns = 200, 150, 100
    x = np.linspace(0, 50, nq)
    y = np.linspace(0, 40, nr)
    z = np.linspace(0, 20, ns)
    
    # Depth-dependent velocity
    Z = z[np.newaxis, np.newaxis, :]
    vs = 1.5 + 0.1 * Z  # km/s, increases with depth
    vp = 1.73 * vs      # Poisson ratio
    
    Qs, Qp = create_q_model(vs, vp, output_prefix='test_q')
    
    print("Created Q model files!")
    print(f"Qs range: [{Qs.min():.1f}, {Qs.max():.1f}]")
    print(f"Qp range: [{Qp.min():.1f}, {Qp.max():.1f}]")
```

---

## 7. File Format Summary

| Aspect | Specification |
|--------|---------------|
| **Format** | Binary, Fortran unformatted |
| **Byte Order** | Native (machine-dependent) |
| **Data Type** | `real*8` (64-bit double precision) |
| **Array Order** | Column-major (k varies fastest) |
| **Header** | None (pure data) |
| **Extensions** | `.bin` or `.dat` |
| **Compression** | None (future: HDF5 with compression) |

### File Size Calculation

```
File size (bytes) = nq × nr × ns × 8 bytes

Example:
  Grid: 500 × 400 × 200 = 40,000,000 points
  Size: 40M × 8 bytes = 320 MB per file
  
  Total for Qs + Qp + gamma: ~960 MB
```

---

## 8. Testing File I/O

### Fortran Test Program

```fortran
program test_q_io
   use common, only : wp
   use withers_io, only : read_q_from_file, write_q_to_file
   
   implicit none
   
   integer, parameter :: nq = 100, nr = 80, ns = 60
   real(kind=wp), dimension(nq,nr,ns) :: Qs_write, Qs_read
   integer :: i, j, k
   
   ! Create test data
   do k = 1, ns
      do j = 1, nr
         do i = 1, nq
            Qs_write(i,j,k) = 50.0_wp + 2.0_wp * real(k, wp)
         end do
      end do
   end do
   
   ! Write to file
   call write_q_to_file('test_Qs.bin', nq, nr, ns, Qs_write)
   
   ! Read back
   call read_q_from_file('test_Qs.bin', nq, nr, ns, Qs_read)
   
   ! Verify
   if (maxval(abs(Qs_write - Qs_read)) < 1e-10_wp) then
      print *, "SUCCESS: Q file I/O works correctly!"
   else
      print *, "FAILURE: Q values don't match!"
   end if
   
end program test_q_io
```

---

## 9. Troubleshooting

### Common Issues

**Issue 1: Wrong array dimensions**
```
Error: Shape mismatch when reading Q file
Solution: Verify file size matches nq×nr×ns×8 bytes
```

**Issue 2: Byte order problems (rare)**
```
Error: Q values are garbage (e.g., 1e308, NaN)
Solution: File created on different architecture
         Convert using Python: array.byteswap().tofile()
```

**Issue 3: File not found**
```
Error: Cannot open Q file
Solution: Check path is relative to run directory
         Use absolute paths if needed
```

**Issue 4: Invalid Q values**
```
Warning: Q values outside valid range [15, 10000]
Solution: Code automatically clamps, but check input data
```

---

## 10. Best Practices

1. **Always validate Q files** after reading
2. **Use consistent byte order** within your workflow
3. **Document Q model provenance** (how it was created)
4. **Test with small models first** (e.g., 10×10×10)
5. **Consider HDF5 for large models** (future implementation)
6. **Keep backup copies** of Q model files
7. **Version control Q models** (use git-lfs for binary files)

---

## Appendix: Complete File Creation Script

See the Python script above in Section 6 for a complete, working example of creating Q model files from scratch or from existing velocity models.
