# rosetta_distribution_opf.jl
 
## Overview

This README describes the two scripts included in this folder:

- `statcom_dclink_scenarios.jl`  
    Explores STATCOM control and device-level behavior under a set of DC‑link and network scenarios. Intended use: study steady‑state interactions of an inverter‑based STATCOM with the distribution network while varying DC‑link parameters, neutral conductor rating, and scenario conditions.

- `SOP_inductionmotors.jl`  
    Computes Soft Open Points (SOPs) and steady operating characteristics for induction motor loads used in distribution studies. Intended use: produce operating points, load models, and parameter sweeps for integration into OPF and dynamic studies.

## Contents (what each script does)

- statcom_dclink_scenarios.jl
    - Defines a set of STATCOM DC‑link ripple power and neutral leg scenarios (2w ripple power, neutral leg rating, input source conditions).
    - Loads a distribution test case and places a STATCOM/inverter model at configured buses.
    - Runs steady‑state optimisation problem for each scenario.
    - Exports results: time series, summary tables, and diagnostic plots (voltage, DC‑link 2w ripple power, currents).
    - Supports parameter sweeps and batch runs for sensitivity analysis.

- SOP_inductionmotors.jl
    - Loads motor nameplate/parameter and derating data، SOP parameters and network operating conditions.
    - Solves for minimising motor derating and DC-link ripple power on the SOP.
    - Outputs linearized models or lookup tables suitable for use in OPF studies.
    - Exports results as CSV/JSON and optionally generates summary plots.

## Usage

1. Install project dependencies:
     - Use Julia 1.6+ (recommended).
     - From the repository root run:
         julia --project=@. -e 'import Pkg; Pkg.instantiate()'

2. Run a script:
     - From the script directory:
         julia --project=@. statcom_dclink_scenarios.jl
         julia --project=@. SOP_inductionmotors.jl

     - Both scripts expose configuration blocks near the top (scenario lists, file paths, solver options). Edit those blocks to change inputs or enable/disable plotting/output.

3. Output locations:
     - Results and figures are written into an `Figures/` subfolder (created on demand).
     - Summary CSV/JSON files and per‑scenario plots are saved with descriptive filenames.

## Inputs and Outputs

- Inputs
    - Network/testcase OpenDSS files.
    - Motor parameter tables (CSV/JSON) for SOP_induction motors.
    - Scenario definitions inside the scripts or external config files if present.

- Outputs
    - Time series (CSV), summary tables (CSV/JSON), and plots (PNG/PDF).

## Dependencies

- Julia 1.6+
- Project uses a Project.toml/Manifest.toml — install with Pkg.instantiate().
- Typical packages used: DifferentialEquations.jl, DataFrames.jl, CSV.jl, Plots.jl (actual list in Project.toml).

## Troubleshooting

- If packages fail to load, run Pkg.instantiate() from the repository root with the correct Julia project.
- Check the configuration blocks at the top of each script for file paths and solver tolerances when results look unexpected.

For more details on network models, control implementations, or to contribute new scenarios, consult the repository-level documentation or open an issue/pull request.