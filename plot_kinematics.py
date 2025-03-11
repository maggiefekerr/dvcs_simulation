# plot_kinematics.py

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
    Boost 'vec' by the 3-velocity 'boost' in units of c.
    Standard pure boost formula.
    """
    bx, by, bz = boost[0], boost[1], boost[2]
    b2 = bx*bx + by*by + bz*bz
    if b2 < 1e-14:
        # No significant boost
        return vec.copy()
    gamma = 1.0 / math.sqrt(1.0 - b2)
    bp = bx*vec[0] + by*vec[1] + bz*vec[2]   # dot product with spatial part
    gamma2 = (gamma - 1.0)/b2

    # Time component after boost
    Eprime = gamma*(vec[3] - bp)
    # Spatial components
    px = vec[0] + gamma2*bp*bx - gamma*bx*vec[3]
    py = vec[1] + gamma2*bp*by - gamma*by*vec[3]
    pz = vec[2] + gamma2*bp*bz - gamma*bz*vec[3]

    return np.array([px, py, pz, Eprime], dtype=float)

def unit_3vec(v):
    """
    Return the unit vector of v's spatial part [v[0], v[1], v[2]].
    If v is nearly zero, return (0,0,0).
    """
    m = mag_3vec(v)
    if m < 1e-14:
        return np.array([0.0, 0.0, 0.0])
    return v[:3]/m


##############################################################################
# 2) Compute phiTrento for the final-state photon
##############################################################################

def compute_phi_trento_photon(e_4vec, ph_4vec, beam_4vec, target_4vec):
    """
    Computes the 'Trento phi' angle for the outgoing photon in DVCS.

    1) Form q = beam - e  (virtual photon 4-vector).
    2) gN_4vec = q + target_4vec  => total gamma*-nucleon system.
    3) Boost e_4vec and ph_4vec (and q if needed) into that gN frame.
    4) Define reference plane with q_gN_unit x e_gN_unit,
       then find angle between that plane and (q_gN_unit x photon_gN_unit).
    5) Return angle in [0,360) degrees.
    """
    # (1) The virtual photon 4-vector
    q_4vec = sub_4vec(beam_4vec, e_4vec)

    # (2) The gamma*-nucleon system total 4vec
    gN_4vec = add_4vec(q_4vec, target_4vec)
    denom = gN_4vec[3]
    if abs(denom) < 1e-14:
        return 0.0
    gN_boost = -gN_4vec[:3]/denom

    # Boost all relevant final vectors into that system
    e_gN  = boost_4vec(e_4vec, gN_boost)
    ph_gN = boost_4vec(ph_4vec, gN_boost)
    q_gN  = boost_4vec(q_4vec, gN_boost)

    # Unit 3-vectors
    q_gN_unit  = unit_3vec(q_gN)
    e_gN_unit  = unit_3vec(e_gN)
    ph_gN_unit = unit_3vec(ph_gN)

    # reference cross product: vT = q_gN_unit x e_gN_unit
    vT = np.cross(q_gN_unit, e_gN_unit)
    mag_vT = np.linalg.norm(vT)
    if mag_vT < 1e-14:
        return 0.0
    vT_unit = vT/mag_vT

    # "other" cross product with the photon: vTH = q_gN_unit x ph_gN_unit
    vTH = np.cross(q_gN_unit, ph_gN_unit)
    mag_vTH = np.linalg.norm(vTH)
    if mag_vTH < 1e-14:
        return 0.0
    vTH_unit = vTH/mag_vTH

    # cosPhi = vT_unit dot vTH_unit
    cosPhi = np.dot(vT_unit, vTH_unit)
    cosPhi = max(-1.0, min(1.0, cosPhi))  # clamp
    phi = math.acos(cosPhi)  # in [0, pi]

    # sign from triple product to get full [0,2π)
    triple = np.dot(np.cross(e_gN_unit, ph_gN_unit), q_gN_unit)
    # or you might prefer a more direct expression:
    # triple = (e_gN[:3] x ph_gN[:3]) . q_gN_unit
    # If triple < 0 => negative sign => phi -> 2π - phi

    if triple < 0.0:
        phi = 2.0*math.pi - phi

    return math.degrees(phi)


##############################################################################
# 3) Main plot function
##############################################################################

def plot_kinematics(lund_file, beam_energy=10.604):
    # Style
    rcParams['font.family'] = 'serif'
    rcParams['font.size'] = 10
    rcParams['axes.labelsize'] = 12
    rcParams['axes.titlesize'] = 14
    rcParams['xtick.labelsize'] = 10
    rcParams['ytick.labelsize'] = 10
    rcParams['legend.fontsize'] = 10
    rcParams['figure.dpi'] = 300
    rcParams['figure.autolayout'] = True

    # Read LUND lines that have 14 columns
    try:
        with open(lund_file, 'r') as f:
            particle_lines = []
            for line in f:
                parts = line.strip().split()
                if len(parts) == 14:
                    particle_lines.append(line)
        data = np.loadtxt(particle_lines)
    except Exception as e:
        print(f"Error loading file: {e}")
        return

    if data.shape[0] == 0:
        print("No valid particle lines found.")
        return

    # data columns: [index, lifetime, type, PID, parent, daughter, px,py,pz, E, mass, vx, vy, vz]
    # We'll separate e-, p, gamma:
    try:
        electrons = data[data[:, 3] == 11]
        protons   = data[data[:, 3] == 2212]
        photons   = data[data[:, 3] == 22]
    except IndexError:
        print("Could not parse e/p/γ from columns.")
        return

    def calc_momentum(part):
        return np.sqrt(part[:, 6]**2 + part[:, 7]**2 + part[:, 8]**2)
    def calc_theta(part):
        p = calc_momentum(part)
        return np.degrees(np.arccos(part[:, 8]/p))

    e_p = calc_momentum(electrons)
    e_theta = calc_theta(electrons)
    p_p = calc_momentum(protons)
    p_theta = calc_theta(protons)
    gamma_p = calc_momentum(photons)
    gamma_theta = calc_theta(photons)

    # Standard DIS variables from electron
    try:
        # E_e = electrons[:, 9]
        # nu = beamE - E_e
        nu = beam_energy - electrons[:, 9]
        Q2 = 4.*beam_energy*electrons[:, 9] * np.sin(np.radians(e_theta)/2.)**2
        y  = nu/beam_energy
        xB = Q2/(2.*0.938272*nu)
        W  = np.sqrt(0.938272**2 + 2.*0.938272*nu - Q2)
        t  = -(protons[:, 6]**2 + protons[:, 7]**2 + (protons[:, 8] - 0.938272)**2)
        photon_phi_lab = np.degrees(np.arctan2(photons[:, 7], photons[:, 6])) % 360
    except Exception as e:
        print(f"Error calculating kinematics: {e}")
        return

    # We want to compute phiTrento for the photon => need 1:1 e- / photon matching
    ne = len(electrons)
    nph = len(photons)
    # also define 4-vectors for the beam and target
    beam_4 = make_4vec(0.0, 0.0, beam_energy, beam_energy)
    target_4 = make_4vec(0.0, 0.0, 0.0, 0.938272)

    nEvents = min(ne, nph)  # simplistic approach: pair up i-th electron with i-th photon
    phi_trento_vals = []
    for i in range(nEvents):
        e_row = electrons[i]
        ph_row = photons[i]
        e_4  = make_4vec(e_row[6], e_row[7], e_row[8], e_row[9])
        ph_4 = make_4vec(ph_row[6], ph_row[7], ph_row[8], ph_row[9])
        phi_tr = compute_phi_trento_photon(e_4, ph_4, beam_4, target_4)
        phi_trento_vals.append(phi_tr)
    phi_trento_vals = np.array(phi_trento_vals)

    # Plotting
    fig, axs = plt.subplots(2, 7, figsize=(28, 8))

    def plot_hist(ax, dat, xlabel, x_range, bins=50):
        counts, bin_edges = np.histogram(dat, bins=bins, range=x_range)
        centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
        ax.errorbar(centers, counts, yerr=np.sqrt(counts), fmt='o',
                    ms=4, lw=1, capsize=2, color='#1f77b4')
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Counts')
        ax.set_xlim(x_range)
        ax.grid(alpha=0.3)

    # First row
    plot_hist(axs[0,0], e_p, r'$e_p$ (GeV)', (0, 12))
    plot_hist(axs[0,1], e_theta, r'$e_\theta$ (deg)', (0, 90))
    plot_hist(axs[0,2], p_p, r'$p_p$ (GeV)', (0, 4))
    plot_hist(axs[0,3], p_theta, r'$p_\theta$ (deg)', (0, 90))
    plot_hist(axs[0,4], gamma_p, r'$\gamma_p$ (GeV)', (0, 10))
    plot_hist(axs[0,5], gamma_theta, r'$\gamma_\theta$ (deg)', (0, 90))
    # The new phiTrento distribution
    plot_hist(axs[0,6], phi_trento_vals, r'$\phi_{\mathrm{Trento}}^\gamma$ (deg)', (0, 360))

    # Second row
    plot_hist(axs[1,0], Q2, r'$Q^2$ (GeV$^2$)', (0, 12))
    plot_hist(axs[1,1], W, r'$W$ (GeV)', (1, 6))
    plot_hist(axs[1,2], xB, r'$x_B$', (0, 1))
    plot_hist(axs[1,3], -t, r'$-t$ (GeV$^2$)', (0, 1))
    plot_hist(axs[1,4], photon_phi_lab, r'$\phi_{\mathrm{lab}}^\gamma$ (deg)', (0, 360))
    axs[1,5].axis('off')
    axs[1,6].axis('off')

    fig.suptitle("Kinematics with Photon Trento $\phi$")
    plt.subplots_adjust(wspace=0.3, hspace=0.4)
    outfile = lund_file.replace('.dat', '_plots.png')
    plt.savefig(outfile)
    plt.close()
    print(f"Successfully created: {outfile}")


def main():
    parser = argparse.ArgumentParser(
        description='Plot kinematics from Lund file, including final-state photon Trento phi',
        add_help=False
    )
    parser.add_argument('input_file', nargs='?', help='Input .dat file')
    parser.add_argument('-b', '--beam-energy', type=float, default=10.604,
                        help='Beam energy in GeV (default: 10.604)')
    
    if len(sys.argv) == 1:
        print("Usage: python plot_kinematics.py INPUT_FILE [OPTIONS]")
        print("  -b, --beam-energy <float>  (default=10.604)")
        return

    args = parser.parse_args()
    if not args.input_file:
        print("Error: no input file specified.")
        return
    
    plot_kinematics(args.input_file, args.beam_energy)


if __name__ == '__main__':
    main()