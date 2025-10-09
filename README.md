# Energy Budget Issues: High-Fidelity Simulations Using Basilisk C

[![Issues](https://img.shields.io/github/issues/comphy-lab/energy-budget-issues)](https://github.com/comphy-lab/energy-budget-issues/issues)
[![License](https://img.shields.io/github/license/comphy-lab/energy-budget-issues)](https://github.com/comphy-lab/energy-budget-issues/blob/main/LICENSE)
[![Last Commit](https://img.shields.io/github/last-commit/comphy-lab/energy-budget-issues)](https://github.com/comphy-lab/energy-budget-issues/commits/main)
[![Basilisk](https://img.shields.io/badge/Basilisk-Compatible-green)](http://basilisk.fr/)
[![Research](https://img.shields.io/badge/Research-Multiphase%20Flow-blue)](https://comphy-lab.org)
![CodeRabbit Pull Request Reviews](https://img.shields.io/coderabbit/prs/github/comphy-lab/energy-budget-issues?utm_source=oss&utm_medium=github&utm_campaign=comphy-lab%2Fenergy-budget-issues&labelColor=171717&color=FF570A&link=https%3A%2F%2Fcoderabbit.ai&label=CodeRabbit+Reviews)

Repository for investigating energy budget issues in multiphase flow simulations using Basilisk C.

## Project Description

This repository focuses on investigating energy budget issues in multiphase flow simulations using Basilisk C. The project examines energy conservation in various fluid dynamics scenarios including drop impact, bubble dynamics, and interface evolution with complex physics such as viscoelasticity and thermal effects.

### Key Research Areas

- **Energy Conservation Analysis** Tracking kinetic, surface, and dissipation energy components
- **Multiphase Flow Dynamics** Drop impact, bubble bursting, and interface evolution
- **Complex Physics Implementation** Viscoelastic fluids, thermal effects, and non-Newtonian rheology
- **High-Fidelity Simulations** Using adaptive mesh refinement and conservative schemes
- **Validation and Benchmarking** Comparison with experimental data and theoretical predictions

### Prerequisites for Contributors

- **Experience with Basilisk C** or similar CFD frameworks
- **Understanding of fluid mechanics** particularly multiphase flows
- **Knowledge of energy conservation principles** in continuum mechanics
- **Programming skills in C** and familiarity with parallel computing

## Repository Structure

```
├── basilisk/src/: Core Basilisk CFD library (reference only)
├── postProcess/: Project-specific post-processing tools
├── src-local/: Custom header files extending Basilisk functionality
├── simulationCases/: Simulation cases with energy budget analysis
```

## Installation

Basilisk requires a C99-compliant compiler and GNU make. Installation can be done in two ways:

### Global Installation

#### Using darcs (recommended)
```bash
sudo apt install darcs make gawk
darcs clone http://basilisk.fr/basilisk
cd basilisk/src

# For Linux/Ubuntu users (preferred mode of operation)
ln -s config.gcc config

# For Mac users
# ln -s config.osx config

make
```

#### Using a tarball
```bash
wget http://basilisk.fr/basilisk/basilisk.tar.gz
tar xzf basilisk.tar.gz
cd basilisk/src

# For Linux/Ubuntu users (preferred mode of operation)
ln -s config.gcc config

# For Mac users
# ln -s config.osx config

make
```

#### Add to your shell configuration (.bashrc or .zshrc)
```bash
echo "export BASILISK=$PWD" >> ~/.bashrc
echo 'export PATH=$PATH:$BASILISK' >> ~/.bashrc
```
Or for zsh users:
```bash
echo "export BASILISK=$PWD" >> ~/.zshrc
echo 'export PATH=$PATH:$BASILISK' >> ~/.zshrc
```

### Repository Level Installation

For project-specific installations, you can use the provided `reset_install_requirements.sh` script which:
- Installs Basilisk within your project directory
- Sets up environment variables locally (in `.project_config`)
- Automatically detects your OS (MacOS or Linux) and uses appropriate configuration
- Verifies the installation

#### Basic usage:
```bash
# Run the script to install or use existing installation
./reset_install_requirements.sh

# For a fresh installation (removes existing one if present)
./reset_install_requirements.sh --hard

# Load the environment settings for your current shell session
source .project_config
```

The script will create a `.project_config` file in your project root with the necessary environment variables. This approach avoids modifying your global shell configuration and keeps the Basilisk installation contained within your project.

### Windows Subsystem for Linux (WSL) Compatibility

Testing on WSL is currently incomplete. In principle, the Linux installation instructions should work for WSL environments. If you encounter any issues while installing or running Basilisk on WSL, please report them by [opening a bug report](https://github.com/comphy-lab/energy-budget-issues/issues/new?template=bug_report.md&labels=bug,wsl).

### Complete Installation Instructions

For more detailed installation instructions, including configuration for different systems, setting up environment variables, installing additional dependencies, and optional libraries, please refer to the official installation guide at [http://basilisk.fr/src/INSTALL](http://basilisk.fr/src/INSTALL).

### Running the codes

To use the make file do:
```bash
CFLAGS=-DDISPLAY=-1 make NAME-of-File.tst
```

## Reporting Issues and Feature Requests

We use GitHub Issues to track bugs, feature requests, and example requests for this course. When creating an issue, please select the appropriate template to help us address your needs efficiently.

### Issue Templates

#### Bug Report:
[Report here](https://github.com/comphy-lab/energy-budget-issues/issues/new?template=bug_report.md)

- For problems with installation, compilation, or running code. 
Please include:
- Detailed description of the issue
- Your environment (OS, compiler version)
- Steps to reproduce
- Expected vs. actual behavior
- Error messages and logs
- Code snippets or files that demonstrate the issue

#### Feature/Topic Request:
[Report here](https://github.com/comphy-lab/energy-budget-issues/issues/new?template=feature_request.md)
- For requesting specific topics or functionality
- Coverage of specific topics in the course
- New examples or tutorials
- Additional functionality in the codebase
- Improvements to existing materials

#### Example Request:
[Report here](https://github.com/comphy-lab/energy-budget-issues/issues/new?template=example_request.md)
- For requesting specific examples that demonstrate:
- Particular Basilisk features
- Solutions to common problems
- Implementation of specific physics or numerical methods

#### General Question:
[Report here](https://github.com/comphy-lab/energy-budget-issues/issues/new?template=general_question.md)
- For any other questions

### How to Create an Issue

1. Go to the ["Issues" tab](https://github.com/comphy-lab/energy-budget-issues/issues) in the GitHub repository
2. Click the ["New Issue"](https://github.com/comphy-lab/energy-budget-issues/issues/new/choose) button
3. Select the appropriate template from the options
4. Fill in the required information according to the template
5. Add relevant labels if available
6. Submit the issue