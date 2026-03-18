# Area-level multi-type model (tau3 version)

This version keeps the separate scaling parameter `tau_3` for the binomial-specific spatial effect in the joint model,
while retaining the current basis-selection rule based on the top 25% largest absolute eigenvalues.

Files:
- `packages.r`: package loader
- `functions.r`: samplers, utilities, and helper functions
- `build_SD_cleaned.R`: build South Dakota tract-level analysis data
- `build_westnorthc.R`: build West North Central tract-level analysis data
- `empirical study.r`: simulation study
- `empirical results.r`: summary tables and maps from saved simulation output
- `region.r`: empirical tract-level analysis and plots
