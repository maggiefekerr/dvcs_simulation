#!/usr/bin/env python3
import sys
import os
import argparse
import numpy as np
import time
import subprocess

# Add current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from km15gen import genOneEvent

# Constants
M = 0.938272081  # Proton mass in GeV/c²

# Bin calculation constants
x1 = 1/(2*M*8.604)
x2 = 1/(5 - M**2)
x3 = (10.604/8.604 - 1)/M*10.604*(1 - np.cos(np.radians(35)))
x4 = (1 - (4 - M**2)/(2*10.604*M)) / (1 + (4 - M**2)/(2*10.604**2*(1 - np.cos(np.radians(35)))))

# Q² bins
y1, y2, y3, y4, y5 = 1.0, 1.456, 2.510, 4.326, 7.671

# xB bins calculations
c0 = y2/(2*M*8.604)
d0 = 1/(1 + (4 - M**2)/y2)
c1 = np.sqrt(y2*y3)/(2*M*8.604)
d1 = 1/(1 + (4 - M**2)/np.sqrt(y2*y3))
c2 = y3/(2*M*8.604)
d2 = 1/(1 + (4 - M**2)/y3)
c3 = np.sqrt(y3*y4)/(2*M*8.604)
d3 = 1/(1 + (4 - M**2)/np.sqrt(y3*y4))
c4 = y4/(2*M*8.604)
d4 = 1/(1 + (4 - M**2)/y4)
c5 = np.sqrt(y4*y5)/(2*M*8.604)
d5 = 1/(1 + (4 - M**2)/np.sqrt(y4*y5))
c6 = y5/(2*M*8.604)
d6 = 1/(1 + (4 - M**2)/y5)

# Bin definitions
newxBbins = [x1, c0, c1, c2, c3, c4, d2]
newQ2bins = [y1, y2, np.sqrt(y2*y3), y3, np.sqrt(y3*y4), y4, np.sqrt(y4*y5)]
newxBbins2 = [x1, c0, c1, c2, c3, c4, c5, d2, d4]
newQ2bins2 = [y1, 1.2, y2, np.sqrt(y2*y3), y3, np.sqrt(y3*y4), y4, np.sqrt(y4*y5), 7]
newtbins = [0.11, 0.15, 0.25, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 1.79]

def main(args):
    """Main function to handle event generation based on specified model"""
    # Set up environment variables for dvcsgen
    os.environ["CLASDVCS_PDF"] = os.path.join(os.path.dirname(__file__), "dependencies", "dvcsgen")
    os.environ["PATH"] = f"{os.environ['CLASDVCS_PDF']}:{os.environ.get('PATH', '')}"
    
    # Parse input parameters
    Ed = args.Ed
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

    # Common parameters dictionary WITH SEED AND BIN
    params = {
        'ymin': args.ymin,
        'ymax': args.ymax,
        'w2min': args.w2min,
        'rad': int(args.radgen),
        'trig': args.trig,
        'filename': args.fname,
        'Ed': Ed,
        'seed': args.seed,  # Critical fix
        'bin': args.bin      # Critical fix
    }

    # Model-specific handling
    if 'km15' in args.model:
        generate_km15_events(xBmin, xBmax, Q2min, Q2max, tmin, tmax, params)
    elif args.model == 'bh':
        generate_bh_events(xBmin, xBmax, Q2min, Q2max, tmin, tmax, params)
    elif args.model == 'vgg':
        generate_vgg_events(xBmin, xBmax, Q2min, Q2max, tmin, tmax, params)
    else:
        print(f"Unknown model: {args.model}")
        sys.exit(1)

def generate_km15_events(xBmin, xBmax, Q2min, Q2max, tmin, tmax, params):
    """Generate events using the KM15 model"""
    start_time = time.time()
    num_events = 0
    output = []
    
    while num_events < params['trig']:
        result = genOneEvent(
            xBmin, xBmax, Q2min, Q2max, tmin, tmax,
            params['ymin'], params['ymax'], params['w2min'], 0,
            rad=params['rad'], Ed=params['Ed'],
            filename=params['filename'], model='km15'
        )
        if result:
            num_events += 1
            output.append(result)
    
    # Write output file
    with open(f"{params['filename']}.dat", "w") as f:
        f.write("".join(output)[:-1])  # Remove trailing newline
    
    print(f"Generated {num_events} KM15 events in {time.time()-start_time:.2f}s")

def generate_bh_events(xBmin, xBmax, Q2min, Q2max, tmin, tmax, params):
    """Generate Bethe-Heitler events using dvcsgen"""
    mode = 0
    if params['rad']:
        mode = 1 if params.get('fringe') else 0
    
    dvcsgen_path = os.path.join(os.path.dirname(__file__), "dependencies", "dvcsgen", "dvcsgen")
    cmd = [
        dvcsgen_path,
        "--docker",
        "--trig", str(params['trig']),
        "--beam", f"{params['Ed']:.3f}",
        "--x", f"{xBmin:.3f}", f"{xBmax:.3f}",
        "--q2", f"{Q2min:.3f}", f"{Q2max:.3f}",
        "--t", f"{tmin:.3f}", f"{tmax:.3f}",
        "--gpd", "101",
        "--y", f"{params['ymin']:.3f}", f"{params['ymax']:.3f}",
        "--w", f"{params['w2min']:.3f}",
        "--fname", params['filename'],
        "--raster", "0.025",
        "--zpos", "-3",
        "--zwidth", "5",
        "--writef", "2",
        "--globalfit",
        "--ycol", "0.0005",
        "--weight",
        "--seed", f"{1000000*mode + 1000*params.get('bin',0) + params['seed']}",
        "--bh", "1"
    ]
    
    if params['rad']:
        cmd += ["--radgen", "--vv2cut", "0.6", "--delta", "0.1", "--radstable"]
    
    print("Executing BH:", " ".join(cmd))
    subprocess.run(cmd, check=True)

def generate_vgg_events(xBmin, xBmax, Q2min, Q2max, tmin, tmax, params):
    """Generate VGG model events using dvcsgen"""
    mode = 3 if params['rad'] else 5
    if params.get('fringe'):
        mode += 1
    
    dvcsgen_path = os.path.join(os.path.dirname(__file__), "dependencies", "dvcsgen", "dvcsgen")
    cmd = [
        dvcsgen_path,
        "--docker",
        "--trig", str(params['trig']),
        "--beam", f"{params['Ed']:.3f}",
        "--x", f"{xBmin:.3f}", f"{xBmax:.3f}",
        "--q2", f"{Q2min:.3f}", f"{Q2max:.3f}",
        "--t", f"{tmin:.3f}", f"{tmax:.3f}",
        "--gpd", "101",
        "--y", f"{params['ymin']:.3f}", f"{params['ymax']:.3f}",
        "--w", f"{params['w2min']:.3f}",
        "--fname", params['filename'],
        "--raster", "0.025",
        "--zpos", "-3",
        "--zwidth", "5",
        "--writef", "2",
        "--globalfit",
        "--ycol", "0.0005",
        "--weight",
        "--seed", f"{1000000*mode + 1000*params.get('bin',0) + params['seed']}",
        "--bh", "3"
    ]
    
    if params['rad']:
        cmd += ["--radgen", "--vv2cut", "0.6", "--delta", "0.1", "--radstable"]
    
    print("Executing VGG:", " ".join(cmd))
    subprocess.run(cmd, check=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="DVCS Event Generator",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Core parameters
    parser.add_argument("-Ed", "--Ed", type=float, default=10.604,
                        help="Beam energy in GeV")
    parser.add_argument("-trig", "--trig", type=int, default=1,
                        help="Number of events to generate")
    parser.add_argument("-fname", "--fname", type=str, default="output",
                        help="Base filename for output files")
    parser.add_argument("-model", "--model", type=str, default='km15',
                        choices=['km15', 'bh', 'vgg'],
                        help="Physics model to use for generation")
    
    # Bin selection
    parser.add_argument("-bin", "--bin", type=int, default=0,
                        help="Use predefined bin scheme (0=disable)")
    parser.add_argument("-fringe", "--fringe", action='store_true',
                        help="Use fringe bin scheme with -bin")
    
    # Kinematic ranges (used when -bin=0)
    parser.add_argument("-xBmin", "--xBmin", type=float, default=0.05,
                        help="Minimum xB when not using bins")
    parser.add_argument("-xBmax", "--xBmax", type=float, default=0.75,
                        help="Maximum xB when not using bins")
    parser.add_argument("-Q2min", "--Q2min", type=float, default=0.9,
                        help="Minimum Q² when not using bins")
    parser.add_argument("-Q2max", "--Q2max", type=float, default=11.0,
                        help="Maximum Q² when not using bins")
    parser.add_argument("-tmin", "--tmin", type=float, default=0.085,
                        help="Minimum |t| when not using bins")
    parser.add_argument("-tmax", "--tmax", type=float, default=1.79,
                        help="Maximum |t| when not using bins")
    
    # Additional cuts
    parser.add_argument("-ymin", "--ymin", type=float, default=0.19,
                        help="Minimum y value")
    parser.add_argument("-ymax", "--ymax", type=float, default=0.85,
                        help="Maximum y value")
    parser.add_argument("-w2min", "--w2min", type=float, default=3.61,
                        help="Minimum W² value")
    
    # Special options
    parser.add_argument("-radgen", "--radgen", action='store_true',
                        help="Enable radiative effects")
    parser.add_argument("-seed", "--seed", type=int, default=0,
                        help="Random seed (0=time-based)")
    
    args = parser.parse_args()
    main(args)