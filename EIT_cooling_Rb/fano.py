import os
import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, steadystate, expect, clebsch

# =============================================================================
# 1. PHYSICAL PARAMETERS & CONFIGURATION
# =============================================================================
gamma = 1.0  # Normalized decay rate
MHz = 1 / 6.067  # Scaling factor to match physical frequencies

# Hyperfine Energies for 87Rb (5P3/2 excited states) relative to F'=2
E_e3 = +266.65 * MHz  # F'=3
E_e2 = 0.0  # F'=2 (Reference)
E_e1 = -156.94 * MHz  # F'=1
E_e0 = -229.16 * MHz  # F'=0

# Trap frequency
nu = 0.2 * gamma

# Laser Parameters (BLUE DETUNED)
Delta_c = +15.0 * gamma  # Positive = Blue Detuned
Delta_p_center = Delta_c  # Probe scans around the coupler detuning

Omega_c_amp = np.sqrt(4 * np.abs(Delta_c) * nu)
Omega_p_amp = 0.15 * gamma

# State Mapping (Total 24 levels)
atom_labels = []
for f, m in [(1, m) for m in range(-1, 2)]:
    atom_labels.append(("g1", f, m))  # 3 states
for f, m in [(2, m) for m in range(-2, 3)]:
    atom_labels.append(("g2", f, m))  # 5 states
for f, m in [(0, m) for m in range(0, 1)]:
    atom_labels.append(("e0", f, m))  # 1 state
for f, m in [(1, m) for m in range(-1, 2)]:
    atom_labels.append(("e1", f, m))  # 3 states
for f, m in [(2, m) for m in range(-2, 3)]:
    atom_labels.append(("e2", f, m))  # 5 states
for f, m in [(3, m) for m in range(-3, 4)]:
    atom_labels.append(("e3", f, m))  # 7 states

N_atom = len(atom_labels)  # 24


# Helper function for Clebsch-Gordan coefficients
def safe_clebsch(j1, j2, j3, m1, m2, m3):
    if not (abs(j1 - j2) <= j3 <= (j1 + j2)):
        return 0.0
    if m1 + m2 != m3:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    return float(clebsch(j1, j2, j3, m1, m2, m3))


# =============================================================================
# 2. SIMULATION & SPECTRUM CALCULATION
# =============================================================================
def calculate_eit_profile():
    print("Starting EIT profile calculation for 24-level Rb87...")

    # Indices of specific levels for non-ideality tracking
    g1_idxs = [i for i, l in enumerate(atom_labels) if l[0] == "g1"]
    g2_idxs = [i for i, l in enumerate(atom_labels) if l[0] == "g2"]
    e_idxs = [i for i, l in enumerate(atom_labels) if l[0].startswith("e")]

    # Unwanted ground states (Optical pumping dark states)
    # The ideal Lambda is g2(F=2, m=-2) -> e2(F'=2, m=-2) <- g1(F=1, m=-1)
    ideal_g_probe = [
        i for i, l in enumerate(atom_labels) if l[1] == 2 and l[2] == -2
    ][0]
    unwanted_g2_idxs = [i for i in g2_idxs if i != ideal_g_probe]
    e3_idxs = [i for i, l in enumerate(atom_labels) if l[1] == 3]

    # Frequency scan range
    scan_width = 4.0 * gamma
    det_list = np.linspace(Delta_c - scan_width, Delta_c + scan_width, 300)

    # Data lists
    absorption_total = []
    pop_unwanted_g2 = []
    pop_g1_total = []
    pop_leakage_e3 = []

    # 1. Constant Hamiltonian (Energies + Coupler + Repumper)
    H_const = 0
    for i in e_idxs:
        label, f_prime, m_prime = atom_labels[i]
        # Energy relative to E_e2
        energy = (globals()[f"E_{label}"] - E_e2) - Delta_c
        H_const += energy * (basis(N_atom, i) * basis(N_atom, i).dag())

    # Coupling Laser (sigma-, q = -1)
    pol_c = -1
    for gi in g2_idxs:
        fg, mg = atom_labels[gi][1], atom_labels[gi][2]
        for ei in e_idxs:
            fe, me = atom_labels[ei][1], atom_labels[ei][2]
            c_coeff = safe_clebsch(fg, 1, fe, mg, pol_c, me)
            if abs(c_coeff) > 1e-5:
                V_c = (
                    (Omega_c_amp * c_coeff / 2.0)
                    * basis(N_atom, ei)
                    * basis(N_atom, gi).dag()
                )
                H_const += V_c + V_c.dag()

    # Repumper Laser (F=1 -> F'=2, pi polarization q=0)
    Omega_repump = 0.5 * gamma
    pol_repump = 0
    for gi in g1_idxs:
        fg, mg = atom_labels[gi][1], atom_labels[gi][2]
        for ei in e_idxs:
            fe, me = atom_labels[ei][1], atom_labels[ei][2]
            if fe == 2:  # Target F'=2
                c_coeff = safe_clebsch(fg, 1, fe, mg, pol_repump, me)
                if abs(c_coeff) > 1e-5:
                    V_rep = (
                        (Omega_repump * c_coeff / 2.0)
                        * basis(N_atom, ei)
                        * basis(N_atom, gi).dag()
                    )
                    H_const += V_rep + V_rep.dag()

    # 2. Spontaneous Emission Collapse Operators
    c_ops = []
    for ei in e_idxs:
        fe, me = atom_labels[ei][1], atom_labels[ei][2]
        for gi in g1_idxs + g2_idxs:
            fg, mg = atom_labels[gi][1], atom_labels[gi][2]
            for q in [-1, 0, 1]:
                c_vac = safe_clebsch(fg, 1, fe, mg, q, me)
                if abs(c_vac) > 1e-5:
                    c_ops.append(
                        np.sqrt(gamma)
                        * c_vac
                        * (basis(N_atom, gi) * basis(N_atom, ei).dag())
                    )

    # 3. Projection Operators for Diagnostics
    P_e_total = sum([basis(N_atom, i) * basis(N_atom, i).dag() for i in e_idxs])
    P_unwanted_g2 = sum(
        [basis(N_atom, i) * basis(N_atom, i).dag() for i in unwanted_g2_idxs]
    )
    P_g1_total = sum(
        [basis(N_atom, i) * basis(N_atom, i).dag() for i in g1_idxs]
    )
    P_e3_leakage = sum(
        [basis(N_atom, i) * basis(N_atom, i).dag() for i in e3_idxs]
    )

    # 4. Probe Scan
    pol_p = 0  # Pi polarization
    for dp in det_list:
        H_probe = 0
        two_photon_det = dp - Delta_c

        for gi in g2_idxs:
            H_probe += -two_photon_det * (
                basis(N_atom, gi) * basis(N_atom, gi).dag()
            )
            fg, mg = atom_labels[gi][1], atom_labels[gi][2]
            for ei in e_idxs:
                fe, me = atom_labels[ei][1], atom_labels[ei][2]
                c_coeff = safe_clebsch(fg, 1, fe, mg, pol_p, me)
                if abs(c_coeff) > 1e-5:
                    V_p = (
                        (Omega_p_amp * c_coeff / 2.0)
                        * basis(N_atom, ei)
                        * basis(N_atom, gi).dag()
                    )
                    H_probe += V_p + V_p.dag()

        # Solve for Steady State
        rho_ss = steadystate(H_const + H_probe, c_ops)

        # Collect data
        absorption_total.append(expect(P_e_total, rho_ss))
        pop_unwanted_g2.append(expect(P_unwanted_g2, rho_ss))
        pop_g1_total.append(expect(P_g1_total, rho_ss))
        pop_leakage_e3.append(expect(P_e3_leakage, rho_ss))

    # =========================================================================
    # 3. PLOTTING & DIAGNOSTICS
    # =========================================================================
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

    # --- Plot 1: Fano Profile & Sidebands ---
    ax1.plot(
        det_list,
        absorption_total,
        "k-",
        lw=2,
        label="Total Probe Absorption (Fano profile)",
    )
    # Highlight EIT Resonance
    ax1.axvline(
        Delta_c,
        color="purple",
        ls="--",
        alpha=0.7,
        label=r"EIT Resonance ($\delta = 0$)",
    )

    # Explicit Sidebands for Blue Detuning
    # Red Sideband: Delta_p = Delta_c - nu (Cools the system)
    # Blue Sideband: Delta_p = Delta_c + nu (Heats the system)
    ax1.axvline(
        Delta_c - nu,
        color="red",
        lw=2,
        ls="-.",
        label=r"Red Sideband ($\Delta_c - \nu$)",
    )
    ax1.axvline(
        Delta_c + nu,
        color="blue",
        lw=2,
        ls="-.",
        label=r"Blue Sideband ($\Delta_c + \nu$)",
    )

    ax1.set_ylabel("Total Excited State Population")
    ax1.set_title("EIT / Fano Absorption Profile for 24-level $^{87}Rb$")
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="upper right")

    # --- Plot 2: Trapping in ground states ---
    ax2.plot(
        det_list,
        pop_g1_total,
        "r--",
        label="Pop. in F=1 (Survived the repumper)",
    )
    ax2.plot(
        det_list,
        pop_unwanted_g2,
        "g-",
        label="Pop. in unwanted $m_F$ of F=2",
    )
    ax2.set_ylabel("Ground State Population")
    ax2.set_title("Non-Ideality 1: Optical Pumping into Spurious Dark States")
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # --- Plot 3: Off-resonant leakage ---
    ax3.plot(
        det_list,
        pop_leakage_e3,
        "m-",
        label="Off-resonant excitation to $F'=3$",
    )
    ax3.set_xlabel(r"Probe Detuning $\Delta_p$ ($\gamma$ units)")
    ax3.set_ylabel("Population")
    ax3.set_title(
        "Non-Ideality 2: Leakage & Light Shifts caused by other $F'$ levels"
    )
    ax3.grid(True, alpha=0.3)
    ax3.legend()

    plt.tight_layout()

    # Save image
    if not os.path.exists("images"):
        os.makedirs("images")
    plt.savefig("images/eit_24level_diagnostics.png", dpi=300)
    print("Plot saved to images/eit_24level_diagnostics.png")
    plt.show()


if __name__ == "__main__":
    calculate_eit_profile()