#!/usr/bin/env python3
import sys
import os
import argparse
import numpy as np
import time
import subprocess
import glob
import math

# Add current directory to Python path (for km15gen)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from km15gen import genOneEvent

# Constants
M = 0.938272081  # Proton mass in GeV/c²
me = 0.000511    # Electron mass in GeV

##############################################################################
# Post-processing: Fix LUND file energy and mass columns for VGG/BH events
##############################################################################
def fix_lund_file(filename):
    """
    Reads the given LUND file and for each particle line (14 columns)
    recalculates column 10 (energy) and column 11 (mass) using:
       E = sqrt(px^2 + py^2 + pz^2 + m^2)
    with m set according to the particle's PDG code:
      - For electrons (PID ±11): use m = 0.000511 GeV
      - For protons (PID 2212): use m = 0.938272081 GeV
      - For photons (PID 22): use m = 0.0 GeV (and energy = |p|)
    All other lines (including header lines with 10 columns) are left unchanged.
    The file is overwritten with the fixed version.
    """
    fixed_lines = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            # Check if this is a particle line (should have 14 columns)
            if len(parts) == 14:
                try:
                    # Extract momentum components (columns 7, 8, 9)
                    px = float(parts[6])
                    py = float(parts[7])
                    pz = float(parts[8])
                    p_mod = math.sqrt(px*px + py*py + pz*pz)
                    # Get PDG code from column 4
                    pid = int(parts[3])
                    if abs(pid) == 11:
                        mass = me
                    elif pid == 2212:
                        mass = M
                    elif pid == 22:
                        mass = 0.0
                    else:
                        # For other particles, leave the values unchanged
                        # (or you can add more cases if desired)
                        mass = float(parts[10])
                    
                    energy = math.sqrt(p_mod**2 + mass**2)
                    # Format energy and mass in exponential format with 7 decimals
                    parts[9] = f"{energy: .7E}"
                    parts[10] = f"{mass: .7E}"
                except Exception as err:
                    print(f"Error processing line: {line}\nError: {err}")
            # Reconstruct the line (join with a space)
            fixed_lines.append(" ".join(parts))
    # Write back the fixed file
    with open(filename, 'w') as f:
        for l in fixed_lines:
            f.write(l + "\n")
    print(f"Post-processed LUND file '{filename}' to fix energy and mass columns.")

##############################################################################
# Model-specific event generators
##############################################################################
def generate_km15_events(params):
    """Generate events using the KM15 model."""
    start_time = time.time()
    output = []
    while len(output) < params['trig']:
        result = genOneEvent(
            params['xBmin'], params['xBmax'],
            params['Q2min'], params['Q2max'],
            params['tmin'], params['tmax'],
            params['ymin'], params['ymax'],
            params['w2min'], 0,
            rad=params['rad'], Ed=params['beam'],
            filename=params['fname'], model='km15'
        )
        if result:
            output.append(result)
    with open(f"{params['fname']}.dat", "w") as f:
        # Remove the final newline if needed
        f.write("".join(output)[:-1])
    print(f"Generated {len(output)} KM15 events in {time.time()-start_time:.2f}s")

def generate_bh_events(params):
    """Generate events using the BH model via dvcsgen."""
    mode = 0
    if params['rad']:
        mode = 1 if params.get('fringe') else 0
    dvcsgen_path = os.path.join(os.path.dirname(__file__), "dependencies", "dvcsgen", "dvcsgen")
    cmd = [
        dvcsgen_path,
        "--trig", str(params['trig']),
        "--beam", f"{params['beam']:.3f}",
        "--x", f"{params['xBmin']:.3f}", f"{params['xBmax']:.3f}",
        "--q2", f"{params['Q2min']:.3f}", f"{params['Q2max']:.3f}",
        "--t", f"{params['tmin']:.3f}", f"{params['tmax']:.3f}",
        "--gpd", "101",
        "--y", f"{params['ymin']:.3f}", f"{params['ymax']:.3f}",
        "--w", f"{params['w2min']:.3f}",
        "--file", params['fname'],
        "--raster", "0.025",
        "--zpos", "-3",
        "--zwidth", "5",
        "--writef", "2",
        "--ycol", "0.0005",
        "--weight",
        "--seed", f"{1000000*mode + 1000*params.get('bin',0) + params['seed']}",
        "--bh", "1"
    ]
    if params['rad']:
        cmd += ["--radgen", "--vv2cut", "0.6", "--delta", "0.1", "--radstable"]
    print("Executing BH:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    # Handle filename truncation (dvcsgen uses the first 8 characters)
    base_prefix = params['fname'][:7]
    generated_files = glob.glob(f"{base_prefix}*")
    if generated_files:
        generated_files.sort(key=os.path.getmtime)
        newest_file = generated_files[-1]
        os.rename(newest_file, f"{params['fname']}.dat")
        print(f"Renamed {newest_file} to {params['fname']}.dat")
    # Fix energy and mass columns in the output file
    fix_lund_file(f"{params['fname']}.dat")

def generate_vgg_events(params):
    """Generate events using the VGG model via dvcsgen."""
    mode = 3 if params['rad'] else 5
    if params.get('fringe'):
        mode += 1
    dvcsgen_path = os.path.join(os.path.dirname(__file__), "dependencies", "dvcsgen", "dvcsgen")
    cmd = [
        dvcsgen_path,
        "--trig", str(params['trig']),
        "--beam", f"{params['beam']:.3f}",
        "--x", f"{params['xBmin']:.3f}", f"{params['xBmax']:.3f}",
        "--q2", f"{params['Q2min']:.3f}", f"{params['Q2max']:.3f}",
        "--t", f"{params['tmin']:.3f}", f"{params['tmax']:.3f}",
        "--gpd", "101",
        "--y", f"{params['ymin']:.3f}", f"{params['ymax']:.3f}",
        "--w", f"{params['w2min']:.3f}",
        "--file", params['fname'],
        "--raster", "0.025",
        "--zpos", "-3",
        "--zwidth", "5",
        "--writef", "2",
        "--ycol", "0.0005",
        "--weight",
        "--seed", f"{1000000*mode + 1000*params.get('bin',0) + params['seed']}",
        "--bh", "3"
    ]
    if params['rad']:
        cmd += ["--radgen", "--vv2cut", "0.6", "--delta", "0.1", "--radstable"]
    print("Executing VGG:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    base_prefix = params['fname'][:7]
    generated_files = glob.glob(f"{base_prefix}*")
    if generated_files:
        generated_files.sort(key=os.path.getmtime)
        newest_file = generated_files[-1]
        os.rename(newest_file, f"{params['fname']}.dat")
        print(f"Renamed {newest_file} to {params['fname']}.dat")
    fix_lund_file(f"{params['fname']}.dat")

##############################################################################
# Main function and argument processing
##############################################################################
def main(args):
    """Main function to parse parameters and select the physics model."""
    # Set environment variables
    os.environ["CLASDVCS_PDF"] = os.path.join(os.path.dirname(__file__), "dependencies", "dvcsgen")
    os.environ["PATH"] = f"{os.environ['CLASDVCS_PDF']}:{os.environ.get('PATH', '')}"

    # Determine kinematic ranges: if a bin is selected, load from the appropriate CSV.
    if args.bin:
        bin_file = "fringe_bin_scheme.csv" if args.fringe else "bin_scheme.csv"
        bin_scheme = np.loadtxt(os.path.join(os.path.dirname(__file__), bin_file), delimiter=',')
        xBmin, xBmax, Q2min, Q2max, tmin, tmax = bin_scheme[args.bin - 1]
    else:
        xBmin = args.xBmin
        xBmax = args.xBmax
        Q2min = args.Q2min
        Q2max = args.Q2max
        tmin = args.tmin
        tmax = args.tmax

    params = {
        'beam': args.beam,
        'xBmin': xBmin,
        'xBmax': xBmax,
        'Q2min': Q2min,
        'Q2max': Q2max,
        'tmin': tmin,
        'tmax': tmax,
        'ymin': args.ymin,
        'ymax': args.ymax,
        'w2min': args.w2min,
        'rad': int(args.radgen),
        'trig': args.trig,
        'fname': args.fname,
        'seed': args.seed,
        'bin': args.bin,
        'fringe': args.fringe
    }

    model_handlers = {
        'km15': generate_km15_events,
        'bh': generate_bh_events,
        'vgg': generate_vgg_events
    }
    if args.model not in model_handlers:
        print(f"Unknown model: {args.model}")
        sys.exit(1)
    model_handlers[args.model](params)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="DVCS Event Generator",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # Core parameters
    parser.add_argument("-beam", "--beam", type=float, default=10.604,
                        help="Beam energy in GeV")
    parser.add_argument("-model", "--model", choices=['km15','bh','vgg'], default='km15',
                        help="Physics model to use")
    parser.add_argument("-trig", "--trig", type=int, default=1,
                        help="Number of events to generate")
    parser.add_argument("-fname", "--fname", default="output",
                        help="Base filename for output")
    # Bin selection
    parser.add_argument("-bin", "--bin", type=int, default=0,
                        help="Use predefined bin scheme (0=disable)")
    parser.add_argument("-fringe", "--fringe", action='store_true',
                        help="Use fringe bin scheme with -bin")
    # Kinematic ranges
    parser.add_argument("-xBmin", "--xBmin", type=float, default=0.05,
                        help="Minimum xB when not using bins")
    parser.add_argument("-xBmax", "--xBmax", type=float, default=0.75,
                        help="Maximum xB when not using bins")
    parser.add_argument("-Q2min", "--Q2min", type=float, default=0.9,
                        help="Minimum Q² [GeV²] when not using bins")
    parser.add_argument("-Q2max", "--Q2max", type=float, default=11.0,
                        help="Maximum Q² [GeV²] when not using bins")
    parser.add_argument("-tmin", "--tmin", type=float, default=0.085,
                        help="Minimum |t| [GeV²] when not using bins")
    parser.add_argument("-tmax", "--tmax", type=float, default=1.79,
                        help="Maximum |t| [GeV²] when not using bins")
    # Additional cuts
    parser.add_argument("-ymin", "--ymin", type=float, default=0.19,
                        help="Minimum y value")
    parser.add_argument("-ymax", "--ymax", type=float, default=0.85,
                        help="Maximum y value")
    parser.add_argument("-w2min", "--w2min", type=float, default=3.61,
                        help="Minimum W² [GeV²]")
    # Special options
    parser.add_argument("-radgen", "--radgen", action='store_true',
                        help="Enable radiative effects")
    parser.add_argument("-seed", "--seed", type=int, default=0,
                        help="Random seed (0=time-based)")
    args = parser.parse_args()
    main(args)