# dvcs_simulation

## Overview

A Python-based event generator for Deeply Virtual Compton Scattering (DVCS) simulations, supporting multiple models:

- **KM15**
- **VGG**
- **BH**: pure Bethe-Heitler process

## Installation

```bash
git clone https://github.com/tbhayward/dvcs_simulation.git
cd dvcs_simulation
chmod +x install.sh
./install.sh
source ~/.bashrc
```

## Command-Line Options

### Core (Required) Parameters

| Option    | Description                           | Default  |
|-----------|---------------------------------------|----------|
| `--beam`  | Beam energy in GeV                    | 10.604   |
| `--model` | Physics model (`km15`, `vgg`, `bh`)   | km15     |
| `--trig`  | Number of events to generate          | 1        |
| `--fname` | Output filename prefix                | output   |

### Kinematic Ranges

| Option     | Description              | Default | Range      |
|------------|--------------------------|---------|------------|
| `--xBmin`  | Minimum $x_B$            | 0.05    | 0.001–0.99 |
| `--xBmax`  | Maximum $x_B$            | 0.75    | 0.001–0.99 |
| `--Q2min`  | Minimum $Q^2$ [GeV²]     | 0.9     | 0.1–15     |
| `--Q2max`  | Maximum $Q^2$ [GeV²]     | 11.0    | 0.1–15     |
| `--tmin`   | Minimum $t$ [GeV²]       | 0.085   | 0.01–2.0   |
| `--tmax`   | Maximum $t$ [GeV²]       | 1.79    | 0.01–2.0   |

### Additional Cuts

| Option      | Description            | Default |
|-------------|------------------------|---------|
| `--ymin`    | Minimum $y$           | 0.19    |
| `--ymax`    | Maximum $y$           | 0.85    |
| `--w2min`   | Minimum $W^2$ [GeV²]  | 3.61    |

### Advanced Options

| Option      | Description                                  |
|-------------|----------------------------------------------|
| `--bin N`   | Use predefined bin scheme (1–6)              |
| `--fringe`  | Use fringe binning with `--bin`              |
| `--radgen`  | Enable radiative effects                     |
| `--seed`    | Set random seed (0 = automatic)              |

## Usage Examples

### Basic KM15 Generation

```bash
python main.py --model km15 --trig 1000 --fname km15_run1
```

### VGG Model with Custom Kinematics

```bash
python main.py --model vgg --beam 10.2 --xBmin 0.1 --xBmax 0.6 --Q2min 1.0 --Q2max 5.0
```

### BH with Radiative Effects and Binning

```bash
python main.py --model bh --bin 3 --radgen --seed 42
```

### High-Statistics Run with Custom t-Range

```bash
python main.py --model km15 --trig 100000 --tmin 0.1 --tmax 1.5 --fname large_run
```

## Output Files

- `<fname>.dat`: Generated events in Lund format.
- Automatic handling of temporary files from **dvcsgen**.
- Consistent naming across all models.

## Dependencies

- Python 3.6+
- Cython
- NumPy
- SciPy
- Gepard (for VGG/BH models)
- CLAS12 environment (via `module load clas12`)

## Troubleshooting

**Q:** Getting `test.1.dat` instead of `test.dat`?  
**A:** The code automatically renames files. Temporary files from **dvcsgen** will be cleaned up.

**Q:** Installation fails with Cython errors?  
**A:** Ensure you have Cython installed:
```bash
pip3 install --user cython
```

**Q:** `dvcsgen` not found?  
**A:** Run `source ~/.bashrc` after installation and verify that the `CLASDVCS_PDF` environment variable is set.