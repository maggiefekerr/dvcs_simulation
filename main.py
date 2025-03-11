#!/usr/bin/env python3
import sys
import os
import argparse
import numpy as np
import time
import subprocess
import glob

# Add current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from km15gen import genOneEvent

# Constants
M = 0.938272081  # Proton mass in GeV/c²

def main(args):
    """Main function to handle event generation"""
    # Set environment variables
    os.environ["CLASDVCS_PDF"] = os.path.join(os.path.dirname(__file__), "dependencies", "dvcsgen")
    os.environ["PATH"] = f"{os.environ['CLASDVCS_PDF']}:{os.environ.get('PATH', '')}"
    
    # Parse parameters
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
        'filename': args.fname,
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

def generate_km15_events(params):
    """Generate events using KM15 model"""
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
            filename=params['filename'], model='km15'
        )
        if result:
            output.append(result)
    
    with open(f"{params['filename']}.dat", "w") as f:
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
        "--trig", str(params['trig']),
        "--beam", f"{params['beam']:.3f}",  # Changed from Ed to beam
        "--x", f"{params['xBmin']:.3f}", f"{params['xBmax']:.3f}",
        "--q2", f"{params['Q2min']:.3f}", f"{params['Q2max']:.3f}",
        "--t", f"{params['tmin']:.3f}", f"{params['tmax']:.3f}",
        "--gpd", "101",
        "--y", f"{params['ymin']:.3f}", f"{params['ymax']:.3f}",
        "--w", f"{params['w2min']:.3f}",
        "--file", params['filename'],
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
    
    # Handle filename truncation (dvcsgen uses first 8 characters)
    base_prefix = params['filename'][:7]
    generated_files = glob.glob(f"{base_prefix}*")
    if generated_files:
        generated_files.sort(key=os.path.getmtime)
        newest_file = generated_files[-1]
        os.rename(newest_file, f"{params['filename']}.dat")
        print(f"Renamed {newest_file} to {params['filename']}.dat")

def generate_vgg_events(params):
    """Generate VGG model events using dvcsgen"""
    mode = 3 if params['rad'] else 5
    if params.get('fringe'):
        mode += 1
    
    dvcsgen_path = os.path.join(os.path.dirname(__file__), "dependencies", "dvcsgen", "dvcsgen")
    cmd = [
        dvcsgen_path,
        "--trig", str(params['trig']),
        "--beam", f"{params['beam']:.3f}",  # Changed from Ed to beam
        "--x", f"{params['xBmin']:.3f}", f"{params['xBmax']:.3f}",
        "--q2", f"{params['Q2min']:.3f}", f"{params['Q2max']:.3f}",
        "--t", f"{params['tmin']:.3f}", f"{params['tmax']:.3f}",
        "--gpd", "101",
        "--y", f"{params['ymin']:.3f}", f"{params['ymax']:.3f}",
        "--w", f"{params['w2min']:.3f}",
        "--file", params['filename'],
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
    
    # Handle filename truncation (dvcsgen uses first 8 characters)
    base_prefix = params['filename'][:7]
    generated_files = glob.glob(f"{base_prefix}*")
    if generated_files:
        generated_files.sort(key=os.path.getmtime)
        newest_file = generated_files[-1]
        os.rename(newest_file, f"{params['filename']}.dat")
        print(f"Renamed {newest_file} to {params['filename']}.dat")

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