#!/usr/bin/env python3
import sys
import os
import argparse
import numpy as np
import time
import subprocess
import glob
import math

# Add current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from km15gen import genOneEvent

# Constants
M = 0.938272081  # Proton mass in GeV/c²

##############################################################################
# Post-processing: Fix LUND file energy and mass columns 
##############################################################################
def fix_lund_file(filename):
    """
    Reads the given LUND file and for each particle line (14 columns)
    recalculates column 10 (energy) and column 11 (mass) using:
       E = sqrt(px^2 + py^2 + pz^2 + m^2)
    with m set according to the particle's PDG code:
      - For electrons (PID ±11): use 0.000511 GeV
      - For protons (PID 2212): use 0.938272081 GeV
      - For photons (PID 22): use 0.0 GeV
    Other particle types remain unchanged.
    The file is overwritten with the fixed version.
    """
    fixed_lines = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            # Check for particle lines (should have 14 columns)
            if len(parts) == 14:
                try:
                    # Extract momentum components from columns 7,8,9 (indices 6,7,8)
                    px = float(parts[6])
                    py = float(parts[7])
                    pz = float(parts[8])
                    p_mod = math.sqrt(px*px + py*py + pz*pz)
                    # Get PID from column 4 (index 3)
                    pid = int(parts[3])
                    if abs(pid) == 11:
                        mass = 0.000511
                    elif pid == 2212:
                        mass = M
                    elif pid == 22:
                        mass = 0.0
                    else:
                        mass = float(parts[10])
                    
                    energy = math.sqrt(p_mod**2 + mass**2)
                    # Update column 10 (index 9) and column 11 (index 10)
                    parts[9] = f"{energy: .7E}"
                    parts[10] = f"{mass: .7E}"
                except Exception as err:
                    print(f"Error processing line: {line}\nError: {err}")
            fixed_lines.append(" ".join(parts))
    with open(filename, 'w') as f:
        for l in fixed_lines:
            f.write(l + "\n")
    print(f"Post-processed LUND file '{filename}' to fix energy and mass columns.")

##############################################################################
# Model-specific event generators
##############################################################################
def generate_km15_events(params):
    """Generate events using KM15 model"""
    start_time = time.time()
    output = []
    while len(output) < params['nentries']:
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
        f.write("".join(output)[:-1])
    print(f"Generated {len(output)} KM15 events in {time.time()-start_time:.2f}s")

def generate_bh_events(params):
    """Generate Bethe-Heitler events using dvcsgen"""
    mode = 0
    if params['rad']:
        mode = 1 if params.get('fringe') else 0
    dvcsgen_path = os.path.join(os.path.dirname(__file__), "dependencies", "dvcsgen", "dvcsgen")
    cmd = [
        dvcsgen_path,
        "--trig", str(params['nentries']),
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
    
    # Combine all output files into one
    base_prefix = params['fname'][:7]
    generated_files = glob.glob(f"{base_prefix}*.dat")
    combined_file = f"{params['fname']}.dat"
    with open(combined_file, "w") as outfile:
        for file in sorted(generated_files):
            with open(file, "r") as infile:
                outfile.write(infile.read())
    for file in generated_files:
        os.remove(file)
    print(f"Combined {len(generated_files)} files into {combined_file}")
    fix_lund_file(combined_file)

def generate_vgg_events(params):
    """Generate VGG model events using dvcsgen"""
    mode = 3 if params['rad'] else 5
    if params.get('fringe'):
        mode += 1
    dvcsgen_path = os.path.join(os.path.dirname(__file__), "dependencies", "dvcsgen", "dvcsgen")
    cmd = [
        dvcsgen_path,
        "--trig", str(params['nentries']),
        "--beam", f"{params['beam']:.3f}",
        "--x", f"{params['xBmin']:.3f}", f"{params['xBmax']:.3f}",
        "--q2", f"{params['Q2min']:.3f}", f"{params['Q2max']:.3f}",
        "--t", f"{params['tmin']:.3f}", f"{params['tmax']:.3f}",
        "--gpd", "101",
        # "--gpd", "3",
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
    
    # Combine all output files into one
    base_prefix = params['fname'][:7]
    generated_files = glob.glob(f"{base_prefix}*.dat")
    combined_file = f"{params['fname']}.dat"
    with open(combined_file, "w") as outfile:
        for file in sorted(generated_files):
            with open(file, "r") as infile:
                outfile.write(infile.read())
    for file in generated_files:
        os.remove(file)
    print(f"Combined {len(generated_files)} files into {combined_file}")
    fix_lund_file(combined_file)

##############################################################################
# Main function and argument parsing
##############################################################################
def main(args):
    """Main function to handle event generation."""
    # Set environment variables
    os.environ["CLASDVCS_PDF"] = os.path.join(os.path.dirname(__file__), "dependencies", "dvcsgen")
    os.environ["PATH"] = f"{os.environ['CLASDVCS_PDF']}:{os.environ.get('PATH', '')}"
    
    # Parse kinematic parameters from bin file or arguments
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
        'nentries': args.nentries,
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
    parser.add_argument("-nentries", "--nentries", type=int, default=1,
                        help="Number of entries to generate")
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
    parser.add_argument("-tmin", "--tmin", type=float, default=0.000,
                        help="Minimum |t| [GeV²] when not using bins")
    parser.add_argument("-tmax", "--tmax", type=float, default=1.00,
                        help="Maximum |t| [GeV²] when not using bins")
    # Additional cuts
    parser.add_argument("-ymin", "--ymin", type=float, default=0.19,
                        help="Minimum y value")
    parser.add_argument("-ymax", "--ymax", type=float, default=0.85,
                        help="Maximum y value")
    parser.add_argument("-w2min", "--w2min", type=float, default=3.61,
                        help="Minimum W² [GeV²]")
    # Special options (radgen not really tested yet)
    parser.add_argument("-radgen", "--radgen", action='store_true',
                        help="Enable radiative effects")
    parser.add_argument("-seed", "--seed", type=int, default=0,
                        help="Random seed (0=time-based)")
    args = parser.parse_args()
    main(args)