# plot_kinematics.py
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

def plot_kinematics(lund_file, beam_energy=10.604):
    # Set up scientific style
    rcParams['font.family'] = 'serif'
    rcParams['font.size'] = 10
    rcParams['axes.labelsize'] = 12
    rcParams['axes.titlesize'] = 14
    rcParams['xtick.labelsize'] = 10
    rcParams['ytick.labelsize'] = 10
    rcParams['legend.fontsize'] = 10
    rcParams['figure.dpi'] = 300
    rcParams['figure.autolayout'] = True

    # Load Lund file data (skip header lines with 10 columns)
    try:
        with open(lund_file, 'r') as f:
            particle_lines = []
            for line in f:
                parts = line.strip().split()
                if len(parts) == 14:  # Select only particle lines
                    particle_lines.append(line)
        data = np.loadtxt(particle_lines)
    except Exception as e:
        print(f"Error loading file: {e}")
        return

    # Extract particles 
    try:
        electrons = data[data[:, 3] == 11]   # particle type 11 = electron
        protons = data[data[:, 3] == 2212]   # particle type 2212 = proton
        photons = data[data[:, 3] == 22]     # particle type 22 = photon
    except IndexError:
        print("Invalid file format - could not extract particles")
        return

    # Calculate kinematic variables
    def calc_momentum(part):
        return np.sqrt(part[:, 6]**2 + part[:, 7]**2 + part[:, 8]**2)
    
    def calc_theta(part):
        return np.degrees(np.arccos(part[:, 8]/calc_momentum(part)))
    
    # Electron variables
    e_p = calc_momentum(electrons)
    e_theta = calc_theta(electrons)
    
    # Proton variables
    p_p = calc_momentum(protons)
    p_theta = calc_theta(protons)
    
    # Photon variables
    gamma_p = calc_momentum(photons)
    gamma_theta = calc_theta(photons)
    
    # Calculate DIS variables
    try:
        nu = beam_energy - electrons[:, 9]
        Q2 = 4 * beam_energy * electrons[:, 9] * np.sin(np.radians(e_theta)/2)**2
        y = nu / beam_energy
        xB = Q2 / (2 * 0.938272 * nu)  # proton mass in GeV
        W = np.sqrt(0.938272**2 + 2 * 0.938272 * nu - Q2)
        t = -(protons[:, 6]**2 + protons[:, 7]**2 + (protons[:, 8] - 0.938272)**2)
        photon_phi = np.degrees(np.arctan2(photons[:, 7], photons[:, 6])) % 360
    except Exception as e:
        print(f"Error calculating kinematics: {e}")
        return

    # Create figure
    fig, axs = plt.subplots(2, 6, figsize=(24, 8))
    
    # First row plots
    def plot_hist(ax, data, xlabel, xrange, bins=50):
        counts, bins = np.histogram(data, bins=bins, range=xrange)
        bin_centers = (bins[1:] + bins[:-1])/2
        ax.errorbar(bin_centers, counts, yerr=np.sqrt(counts), 
                   fmt='o', ms=4, lw=1, capsize=2, color='#1f77b4')
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Counts')
        ax.set_xlim(xrange)
        ax.grid(alpha=0.3)
    
    # First row
    plot_hist(axs[0,0], e_p, r'$e_p$ (GeV)', (0, 12))
    plot_hist(axs[0,1], e_theta, r'$e_\theta$ (deg)', (0, 90))
    plot_hist(axs[0,2], p_p, r'$p_p$ (GeV)', (0, 4))
    plot_hist(axs[0,3], p_theta, r'$p_\theta$ (deg)', (0, 90))
    plot_hist(axs[0,4], gamma_p, r'$\gamma_p$ (GeV)', (0, 10))
    plot_hist(axs[0,5], gamma_theta, r'$\gamma_\theta$ (deg)', (0, 90))
    
    # Second row
    plot_hist(axs[1,0], Q2, r'$Q^2$ (GeV$^2$)', (0, 12))
    plot_hist(axs[1,1], W, r'$W$ (GeV)', (1, 6))
    plot_hist(axs[1,2], xB, r'$x_B$', (0, 1))
    plot_hist(axs[1,3], -t, r'$-t$ (GeV$^2$)', (0, 1))
    plot_hist(axs[1,4], photon_phi, r'$\phi$ (deg)', (0, 360))
    axs[1,5].axis('off')
    
    # Adjust spacing and save
    plt.subplots_adjust(wspace=0.3, hspace=0.4)
    output_file = lund_file.replace('.dat', '_plots.png')
    plt.savefig(output_file)
    plt.close()
    print(f"Successfully created: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Plot kinematics from Lund format DVCS simulation data',
        add_help=False
    )
    parser.add_argument('input_file', nargs='?', help='Input .dat file')
    parser.add_argument('-b', '--beam-energy', type=float, default=10.604,
                       help='Beam energy in GeV (default: 10.604)')
    
    if len(sys.argv) == 1:
        print("Usage: python plot_kinematics.py INPUT_FILE [OPTIONS]")
        print("Options:")
        print("  -b, --beam-energy FLOAT   Set beam energy in GeV (default: 10.604)")
        print("  -h, --help                Show this help message")
        return
    
    args = parser.parse_args()
    
    if not args.input_file:
        print("Error: No input file specified")
        print("Usage: python plot_kinematics.py INPUT_FILE [OPTIONS]")
        return
    
    plot_kinematics(args.input_file, args.beam_energy)

if __name__ == '__main__':
    main()