from km15gen import genOneEvent
import argparse
import numpy as np
import time
import os
import subprocess

M = 0.938272081
x1 = 1/2/M/8.604
x2 = 1/(5-M**2)
x3 = (10.604/8.604-1)/M*10.604*(1-np.cos(np.radians(35)))
x4 = (1-(4-M**2)/2/10.604/M/(1+(4-M**2)/2/10.604**2/(1-np.cos(np.radians(35)))))

y1 = 1
y2 = 1.456
y3 = 2.510
y4 = 4.326
y5 = 7.671

c0 = y2/2/M/8.604
d0 = 1/(1+(4-M*M)/y2)
c1 = np.sqrt(y2*y3)/2/M/8.604
d1 = 1/(1+(4-M*M)/np.sqrt(y2*y3))
c2 = y3/2/M/8.604
d2 = 1/(1+(4-M*M)/y3)
c3 = np.sqrt(y3*y4)/2/M/8.604
d3 = 1/(1+(4-M*M)/np.sqrt(y3*y4))
c4 = y4/2/M/8.604
d4 = 1/(1+(4-M*M)/y4)
c5 = np.sqrt(y4*y5)/2/M/8.604
d5 = 1/(1+(4-M*M)/np.sqrt(y4*y5))
c6 = y5/2/M/8.604
d6 = 1/(1+(4-M*M)/y5)

newxBbins = [x1, c0, c1, c2, c3, c4, d2]
newQ2bins = [y1, y2, np.sqrt(y2*y3), y3, np.sqrt(y3*y4), y4, np.sqrt(y4*y5)]
newxBbins2 = [x1, c0, c1, c2, c3, c4, c5, d2, d4]
newQ2bins2 = [y1, 1.2, y2, np.sqrt(y2*y3), y3, np.sqrt(y3*y4), y4, np.sqrt(y4*y5), 7]
newtbins = [0.11, 0.15, 0.25, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 1.79]

def main(args):
    Ed = args.Ed
    if args.bin:
        if args.fringe:
            bin_scheme = np.loadtxt(os.path.join(os.path.dirname(__file__), "fringe_bin_scheme.csv"), delimiter=',')
        else:
            bin_scheme = np.loadtxt(os.path.join(os.path.dirname(__file__), "bin_scheme.csv"), delimiter=',')
        xBmin, xBmax, Q2min, Q2max, tmin, tmax = bin_scheme[args.bin - 1]
    else:
        xBmin = args.xBmin
        xBmax = args.xBmax
        Q2min = args.Q2min
        Q2max = args.Q2max
        tmin = args.tmin
        tmax = args.tmax
    
    ymin = args.ymin
    ymax = args.ymax
    w2min = args.w2min
    rad = int(args.radgen)
    trig = args.trig
    filename = args.fname

    if 'km15' in args.model:
        now = time.time()
        num = 0
        result = ""
        while num < trig:
            this_result = genOneEvent(xBmin, xBmax, Q2min, Q2max, tmin, tmax, ymin, ymax, w2min, 0, rad=rad, Ed=Ed, filename=filename, model=args.model)
            if this_result:
                num += 1
                result += this_result

        later = time.time()
        print(f"The time spent in generating events: {later - now:.3f} s")
        with open(f"{filename}.dat", "w") as fp:
            fp.write(result[:-1])

    elif args.model == 'bh':
        if args.radgen:
            mode = 0
            if args.fringe:
                mode = 1
        else:
            mode = 2
        
        dvcsgen_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dependencies", "dvcsgen", "dvcsgen")
        dvcsgen_commands = [
            dvcsgen_path, "--docker", "--trig", f"{trig}",
            "--beam", f"{Ed:.3f}",
            "--x", f"{xBmin:.3f}", f"{xBmax:.3f}",
            "--q2", f"{Q2min:.3f}", f"{Q2max:.3f}",
            "--t", f"{tmin:.3f}", f"{tmax:.3f}",
            "--gpd", "101", "--y", f"{ymin:.3f}", f"{ymax:.3f}",
            "--w", f"{w2min:.3f}", "--raster", "0.025",
            "--zpos", "-3", "--zwidth", "5", "--writef", "2",
            "--globalfit", "--ycol", "0.0005", "--weight",
            "--seed", f"{1000000*mode + 1000*args.bin + args.seed}"
        ]
        if rad:
            dvcsgen_commands.extend(["--radgen", "--vv2cut", "0.6", "--delta", "0.1", "--radstable"])
        dvcsgen_commands.extend(["--bh", "1"])
        subprocess.run(dvcsgen_commands)

    elif args.model == 'vgg':
        if args.radgen:
            mode = 3
            if args.fringe:
                mode = 4
        else:
            mode = 5
        
        dvcsgen_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dependencies", "dvcsgen", "dvcsgen")
        dvcsgen_commands = [
            dvcsgen_path, "--docker", "--trig", f"{trig}",
            "--beam", f"{Ed:.3f}",
            "--x", f"{xBmin:.3f}", f"{xBmax:.3f}",
            "--q2", f"{Q2min:.3f}", f"{Q2max:.3f}",
            "--t", f"{tmin:.3f}", f"{tmax:.3f}",
            "--gpd", "101", "--y", f"{ymin:.3f}", f"{ymax:.3f}",
            "--w", f"{w2min:.3f}", "--raster", "0.025",
            "--zpos", "-3", "--zwidth", "5", "--writef", "2",
            "--globalfit", "--ycol", "0.0005", "--weight",
            "--seed", f"{1000000*mode + 1000*args.bin + args.seed}"
        ]
        if rad:
            dvcsgen_commands.extend(["--radgen", "--vv2cut", "0.6", "--delta", "0.1", "--radstable"])
        dvcsgen_commands.extend(["--bh", "3"])
        print("Executing:", " ".join(dvcsgen_commands))
        subprocess.run(dvcsgen_commands)

    elif args.model == 'pi0':
        print("Pi0 generator not implemented in this version")
        # Remove or modify pi0 section as needed

    else:
        print(f"Unknown model: {args.model}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="DVCS Simulation Generator", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-Ed", "--Ed", type=float, default=10.604)
    parser.add_argument("-trig", "--trig", type=int, default=1)
    parser.add_argument("-fname", "--fname", type=str, default="km15gen")
    parser.add_argument("-bin", "--bin", type=int, default=0)
    parser.add_argument("-model", "--model", type=str, default='km15')
    parser.add_argument("-xBmin", "--xBmin", type=float, default=0.05)
    parser.add_argument("-xBmax", "--xBmax", type=float, default=0.75)
    parser.add_argument("-Q2min", "--Q2min", type=float, default=0.9)
    parser.add_argument("-Q2max", "--Q2max", type=float, default=11)
    parser.add_argument("-tmin", "--tmin", type=float, default=0.085)
    parser.add_argument("-tmax", "--tmax", type=float, default=1.79)
    parser.add_argument("-ymin", "--ymin", type=float, default=0.19)
    parser.add_argument("-ymax", "--ymax", type=float, default=0.85)
    parser.add_argument("-w2min", "--w2min", type=float, default=3.61)
    parser.add_argument("-radgen", "--radgen", action='store_true')
    parser.add_argument("-fringe", "--fringe", action='store_true')
    parser.add_argument("-seed", "--seed", type=int, default=0)
    args = parser.parse_args()
    main(args)