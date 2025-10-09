# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Basilisk C computational fluid dynamics (CFD) project focused on multiphase flow simulations with complex physics including viscoelasticity, thermal effects, and interface dynamics.

## Build and Development Commands

### Initial Setup
```bash
# Install Basilisk locally in the project directory
./reset_install_requirements.sh

# For a clean reinstallation
./reset_install_requirements.sh --hard

# Load environment variables (required after each new terminal session)
source .project_config
```

### Building Simulations
```bash
# Build and run a test case from simulationCases/
cd simulationCases
CFLAGS=-DDISPLAY=-1 make 3-BurstingBubbles.tst

# Direct compilation with qcc (Basilisk's compiler)
qcc -I../src-local -O2 -Wall -disable-dimensions simulation.c -o simulation -lm

# Run using the provided script
./runCases.sh 3-DropImpactOnPool
```

### Clean Build
```bash
cd simulationCases
./cleanup.sh
```

## Architecture and Code Organization

### Core Components

1. **Custom Physics Headers** (`src-local/`):
   - none so far

2. **Simulation Structure**:
   - Each simulation is a standalone C file that includes relevant Basilisk headers
   - Uses preprocessor directives extensively for configuration
   - Output files are binary format, post-processed separately
   - Volume-of-Fluid (VOF) method for interface tracking

3. **Key Basilisk Patterns**:
   ```c
   // Typical simulation structure
   #include "navier-stokes/centered.h"
   #include "two-phase.h"
   #include "src-local/log-conform-viscoelastic.h"
   
   // Define physical parameters
   double Re = 100;  // Reynolds number
   
   int main() {
     // Initialize grid
     init_grid(1 << LEVEL);
     
     // Set boundary conditions
     // Run simulation
     run();
   }
   ```

### Important Technical Details

- **Compiler**: Uses `qcc` (Basilisk's preprocessor/compiler) which extends C99
- **Parallelization**: Supports MPI for parallel simulations
- **Adaptive Mesh Refinement**: Built-in AMR capabilities
- **Dimension Independence**: Code can run in 2D/3D with minimal changes
- **Conservative Schemes**: Momentum-conserving advection for multiphase flows

### Testing Approach

- Test files have `.tst` extension and are run through Basilisk's test framework
- The Makefile automatically handles test compilation and execution
- Post-processing validation typically done with Python scripts in `postProcess/`

## Development Notes

- Always source `.project_config` before working with Basilisk commands
- Include paths must point to `src-local` for custom headers
- Use `-DDISPLAY=-1` flag for batch/headless runs
- Binary output files can be large; use `cleanup.sh` regularly
- Simulation parameters are typically defined as preprocessor macros or global variables