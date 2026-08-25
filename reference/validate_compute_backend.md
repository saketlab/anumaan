# Validate a selectable Stan compute backend configuration

Validates a compute configuration for CPU or OpenCL execution. This
checks only package/toolchain availability and the structure of the
requested configuration; it does not inspect machine-specific GPU
inventories.

## Usage

``` r
validate_compute_backend(compute = list())
```

## Arguments

- compute:

  Named list. Supported fields are `backend`, `opencl_platform_id`,
  `opencl_device_id`, and `allow_cpu_fallback`.

## Value

A canonical validated list with fields `backend`, `opencl_platform_id`,
`opencl_device_id`, and `allow_cpu_fallback`.
