#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import rcParams

##############################################################################
# 1) Simple 4-vector utilities
##############################################################################
def make_4vec(px, py, pz, E):
    """Return a NumPy array [px, py, pz, E]."""
    return np.array([px, py, pz, E], dtype=float)

def add_4vec(a, b):
    """4-vector addition: a + b."""
    return a + b

def sub_4vec(a, b):
    """4-vector subtraction: a - b."""
    return a - b

def mag_3vec(v):
    """3-vector magnitude of the spatial part (v[0:3])."""
    return np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)

def boost_4vec(vec, boost):
    """
    Boost 'vec' by the 3-velocity 'boost'.
    Standard pure boost formula.
    """
    bx, by, bz = boost[0], boost[1], boost[2]
    b2 = bx*bx + by*by + bz*bz
    if b2 < 1e-14:
        return vec.copy()
    gamma = 1.0 / math.sqrt(1.0 - b2)
    bp = bx*vec[0] + by*vec[1] + bz*vec[2]
    gamma2 = (gamma - 1.0)/b2
    Eprime = gamma*(vec[3] - bp)
    px = vec[0] + gamma2*bp*bx - gamma*bx*vec[3]
    py = vec[1] + gamma2*bp*by - gamma*by*vec[3]
    pz = vec[2] + gamma2*bp*bz - gamma*bz*vec[3]
    return np.array([px, py, pz, Eprime], dtype=float)

def unit_3vec(v):
    """
    Return the unit vector of v's spatial part [v[0], v[1], v[2]].
    If v is nearly zero, return [0,0,0].
    """
    m = mag_3vec(v)
    if m < 1e-14:
        return np.array([0.0, 0.0, 0.0])
    return v[:3] / m

##############################################################################
# 2) Compute phiTrento for the final-state photon
##############################################################################
def compute_phi_trento_photon(e_4vec, ph_4vec, beam_4vec, target_4vec):
    """
    Computes the 'Trento phi' angle for the outgoing photon in DVCS.
    
    1) Form q = beam - e (virtual photon 4-vector).
    2) gN = q + target (total gamma*-nucleon system).
    3) Boost e and photon 4-vectors into the gN frame.
    4) Define the reference plane: vT = (q_unit x e_unit), and compute the angle
       between vT and (q_unit x ph_unit). Use the sign from the triple product.
    5) Return the angle (in degrees) in [0,360).
    """
    q_4vec = sub_4vec(beam_4vec, e_4vec)
    gN_4vec = add_4vec(q_4vec, target_4vec)
    denom = gN_4vec[3]
    if abs(denom) < 1e-14:
        return 0.0
    gN_boost = -gN_4vec[:3] / denom
    
    e_gN = boost_4vec(e_4vec, gN_boost)
    ph_gN = boost_4vec(ph_4vec, gN_boost)
    q_gN = boost_4vec(q_4vec, gN_boost)
    
    q_unit = unit_3vec(q_gN)
    e_unit = unit_3vec(e_gN)
    ph_unit = unit_3vec(ph_gN)
    
    vT = np.cross(q_unit, e_unit)
    mag_vT = np.linalg.norm(vT)
    if mag_vT < 1e-14:
        return 0.0
    vT_unit = vT / mag_vT
    
    vTH = np.cross(q_unit, ph_unit)
    mag_vTH = np.linalg.norm(vTH)
    if mag_vTH < 1e-14:
        return 0.0
    vTH_unit = vTH / mag_vTH
    
    cosPhi = np.dot(vT_unit, vTH_unit)
    cosPhi = max(-1.0, min(1.0, cosPhi))
    phi = math.acos(cosPhi)
    
    triple = np.dot(np.cross(e_unit, ph_unit), q_unit)
    if triple < 0.0:
        phi = 2.0 * math.pi - phi
        
    return math.degrees(phi)

##############################################################################
# 3) Main plot function (supports 1-3 input files with legend options)
##############################################################################
def plot_kinematics(input_files, beam_energy=10.604, legend_labels=None):
    """
    Plots kinematics from 1 to 3 LUND .dat files.
    Top row (6 subplots): e_p, e_theta, p_p, p_theta, gamma_p, gamma_theta.
    Bottom row (6 subplots): y, Q^2, W, x_B, -t, phiTrento (photon).
    
    If multiple datasets are provided, each is plotted with a distinct color:
      - 1 dataset: black
      - 2 datasets: red, blue
      - 3 datasets: red, blue, green
    A legend is added in the top-right of each subplot.
    
    Parameters:
      input_files: list of .dat filenames (1 to 3).
      beam_energy: beam energy in GeV.
      legend_labels: Optional list of legend labels (if not provided, defaults are used).
    """
    # Set style
    rcParams['font.family'] = 'serif'
    rcParams['font.size'] = 10
    rcParams['axes.labelsize'] = 12
    rcParams['axes.titlesize'] = 14
    rcParams['xtick.labelsize'] = 10
    rcParams['ytick.labelsize'] = 10
    rcParams['legend.fontsize'] = 10
    rcParams['figure.dpi'] = 300
    rcParams['figure.autolayout'] = True

    # Define color mapping based on the number of datasets
    color_map = {1: ['black'], 2: ['red', 'blue'], 3: ['red', 'blue', 'green']}
    nfiles = len(input_files)
    if nfiles not in [1, 2, 3]:
        print("Error: Please provide 1 to 3 input files.")
        return
    colors = color_map[nfiles]
    
    # Use default legend labels if not provided
    if legend_labels is None or len(legend_labels) != nfiles:
        legend_labels = [f"Dataset {i+1}" for i in range(nfiles)]
    
    # Initialize lists for data from each file
    all_data = []
    all_electrons = []
    all_protons = []
    all_photons = []
    for fname in input_files:
        try:
            with open(fname, 'r') as f:
                particle_lines = [line for line in f if len(line.strip().split()) == 14]
            data = np.loadtxt(particle_lines)
            all_data.append(data)
            all_electrons.append(data[data[:, 3] == 11])
            all_protons.append(data[data[:, 3] == 2212])
            all_photons.append(data[data[:, 3] == 22])
        except Exception as e:
            print(f"Error processing file {fname}: {e}")
            return

    # Compute phiTrento for each dataset
    phi_trento_all = []
    for i in range(nfiles):
        electrons = all_electrons[i]
        photons = all_photons[i]
        nEvents = min(len(electrons), len(photons))
        beam_4 = make_4vec(0.0, 0.0, beam_energy, beam_energy)
        target_4 = make_4vec(0.0, 0.0, 0.0, 0.938272)
        phi_vals = []
        for j in range(nEvents):
            e_4 = make_4vec(electrons[j][6], electrons[j][7], electrons[j][8], electrons[j][9])
            ph_4 = make_4vec(photons[j][6], photons[j][7], photons[j][8], photons[j][9])
            phi_vals.append(compute_phi_trento_photon(e_4, ph_4, beam_4, target_4))
        phi_trento_all.append(np.array(phi_vals))

    # For each dataset, compute kinematic quantities
    e_p_all = []
    e_theta_all = []
    p_p_all = []
    p_theta_all = []
    gamma_p_all = []
    gamma_theta_all = []
    y_all = []
    Q2_all = []
    W_all = []
    xB_all = []
    t_all = []
    
    for data in all_data:
        def calc_momentum(part):
            return np.sqrt(part[:, 6]**2 + part[:, 7]**2 + part[:, 8]**2)
        def calc_theta(part):
            p = calc_momentum(part)
            return np.degrees(np.arccos(part[:, 8] / p))
        
        electrons = data[data[:, 3] == 11]
        protons   = data[data[:, 3] == 2212]
        photons   = data[data[:, 3] == 22]
        
        e_p_all.append(calc_momentum(electrons))
        e_theta_all.append(calc_theta(electrons))
        p_p_all.append(calc_momentum(protons))
        p_theta_all.append(calc_theta(protons))
        gamma_p_all.append(calc_momentum(photons))
        gamma_theta_all.append(calc_theta(photons))
        
        # DIS kinematics (using electron energy from col 10)
        nu = beam_energy - electrons[:, 9]
        Q2_val = 4.0 * beam_energy * electrons[:, 9] * np.sin(np.radians(calc_theta(electrons))/2.0)**2
        y_val = nu / beam_energy
        xB_val = Q2_val / (2.0 * 0.938272 * nu)
        W_val = np.sqrt(0.938272**2 + 2.0 * 0.938272 * nu - Q2_val)
        t_val = -(protons[:, 6]**2 + protons[:, 7]**2 + (protons[:, 8] - 0.938272)**2)
        
        y_all.append(y_val)
        Q2_all.append(Q2_val)
        W_all.append(W_val)
        xB_all.append(xB_val)
        t_all.append(t_val)
    
    # Create figure with 2 rows x 6 columns
    fig, axs = plt.subplots(2, 6, figsize=(28, 10))
    
    def plot_multiple(ax, data_list, xlabel, x_range, bins=50):
        for d, col, label in zip(data_list, colors, legend_labels):
            counts, bin_edges = np.histogram(d, bins=bins, range=x_range)
            centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
            ax.errorbar(centers, counts, yerr=np.sqrt(counts), fmt='o',
                        ms=4, lw=1, capsize=2, color=col, label=label)
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Counts')
        ax.set_xlim(x_range)
        ax.grid(alpha=0.3)
        ax.legend(loc='upper right')
    
    # Top row plots: e_p, e_theta, p_p, p_theta, gamma_p, gamma_theta
    plot_multiple(axs[0,0], e_p_all, r'$e_p$ (GeV)', (0, 12))
    plot_multiple(axs[0,1], e_theta_all, r'$e_\theta$ (deg)', (0, 90))
    plot_multiple(axs[0,2], p_p_all, r'$p_p$ (GeV)', (0, 4))
    plot_multiple(axs[0,3], p_theta_all, r'$p_\theta$ (deg)', (0, 90))
    plot_multiple(axs[0,4], gamma_p_all, r'$\gamma_p$ (GeV)', (0, 10))
    plot_multiple(axs[0,5], gamma_theta_all, r'$\gamma_\theta$ (deg)', (0, 90))
    
    # Bottom row plots: y, Q2, W, xB, -t, phiTrento
    plot_multiple(axs[1,0], y_all, r'$y$', (0, 1))
    plot_multiple(axs[1,1], Q2_all, r'$Q^2$ (GeV$^2$)', (0, 12))
    plot_multiple(axs[1,2], W_all, r'$W$ (GeV)', (1, 6))
    plot_multiple(axs[1,3], xB_all, r'$x_B$', (0, 1))
    plot_multiple(axs[1,4], t_all, r'$-t$ (GeV$^2$)', (0, 1))
    plot_multiple(axs[1,5], phi_trento_all, r'$\phi$ (deg)', (0, 360))
    
    plt.subplots_adjust(wspace=0.3, hspace=0.4)
    outfile = input_files[0].replace('.dat', '_plots.png')
    plt.savefig(outfile)
    plt.close()
    print(f"Successfully created: {outfile}")

##############################################################################
# Main function and argument parsing
##############################################################################
def main():
    parser = argparse.ArgumentParser(
        description='Plot kinematics from up to three LUND .dat files with legend options',
        add_help=False
    )
    parser.add_argument('input_files', nargs='+', help='Input .dat file(s) (up to 3)')
    parser.add_argument('-b', '--beam-energy', type=float, default=10.604,
                        help='Beam energy in GeV (default: 10.604)')
    parser.add_argument('-l', '--labels', nargs='+',
                        help='Legend labels for each dataset (up to 3)')
    if len(sys.argv) == 1:
        print("Usage: python plot_kinematics.py FILE1.dat [FILE2.dat FILE3.dat] [OPTIONS]")
        print("  -b, --beam-energy <float>  (default=10.604)")
        print("  -l, --labels <label1> [label2 label3]")
        return
    args = parser.parse_args()
    plot_kinematics(args.input_files, args.beam_energy, args.labels)

if __name__ == '__main__':
    main()