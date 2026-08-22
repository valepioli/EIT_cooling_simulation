import numpy as np
from sympy import Rational, sqrt as sqrt_sym, N
from sympy.physics.wigner import wigner_6j, wigner_3j

RUN_NAME = "eit_axial"

# ==========================================
# --- PHYSICAL CONSTANTS & SCALING ---
# ==========================================
gamma = 1.0                # Normalized decay rate
MHz = 1 / 6.067             # Scaling factor (1 MHz in units of Gamma), gamma/2pi = 6.067 MHz

E_e3 = +266.65 * MHz
E_e2 = 0.0
E_e1 = -156.94 * MHz
E_e0 = -229.16 * MHz

# ==========================================
# --- AXIAL SIMULATION PARAMETERS (Conveyor Belt) ---
# ==========================================
# FIX: previously "nu = 0.07 * gamma" -> 424.7 kHz (a ~6% mismatch with the
# nu_z,tip = 450 kHz analytically derived in Sec. "Simulation with Axial
# Experimental Parameters" / Conveyor_belt.pdf). Set nu so that it converts
# to exactly 450 kHz given gamma/2pi = 6.067 MHz.
nu = 0.450 / 6.067 * gamma   # -> 450 kHz exactly
eta = 0.09                  # Axial Lamb-Dicke parameter (probe beam, strong confinement)
n_phon = 2                  # Initial phonon population derived from Tz = 52.8 uK
N_vib = 10                  # Manageable Fock space ceiling for n=2

# ==========================================
# --- MAGNETIC FIELD PARAMETERS ---
# ==========================================
B_field = 4.0              # Gauss
mu_B = 1.399 * MHz         # Bohr magneton in MHz/Gauss
g_g1, g_g2 = -0.5, 0.5

# FIX: corrected scheme (matches Fig. 1 of the thesis / repumper-pi scheme):
#   - probe  (sigma-): |F=2, m=-1> -> |e>   (was pi on |F=2, m=-2> -> |e>)
#   - coupling (sigma-): |F=1, m=-1> -> |e> (unchanged)
#   - repump (pi):     |F=2, m=-2> -> |e>   (was sigma- on |F=2, m=-1> -> |e>)
# Both EIT ground states of the closed Lambda system now sit at m=-1.
m_g1, m_g2 = -1, -1

# ==========================================
# --- LASER PARAMETERS ---
# ==========================================
Delta_c = +13 * gamma
zeeman_offset = (g_g2 * m_g2 - g_g1 * m_g1) * mu_B * B_field
Delta_p_center = Delta_c - zeeman_offset

# FIX: use the *exact* EIT cooling condition, Omega_c^2 = 4*nu*(nu+Delta_c)
# (thesis Eq. "exact_condition_numerical"), instead of only ever using the
# nu << Delta_c approximation Omega_c^2 ~= 4*Delta_c*nu. The exact form
# reduces smoothly to the approximate one in this regime, so it is a strict
# improvement with no other side effects.
Omega_c_amp = np.sqrt(4 * nu * (nu + np.abs(Delta_c)))
Omega_p_amp = 0.09 * gamma
Omega_r_amp = 0.3 * gamma

# FIX: corrected polarizations.
#   pol_p = -1 (sigma-)  probe:    was 0 (pi)
#   pol_c = -1 (sigma-)  coupling: unchanged
#   pol_r =  0 (pi)      repump:   was -1 (sigma-)
pol_p = -1
pol_c = -1
pol_r = 0

# ==========================================
# --- STATE MAPPING ---
# ==========================================
atom_labels = []
for f, m in [(1, m) for m in range(-1, 2)]: atom_labels.append(("g1", f, m))
for f, m in [(2, m) for m in range(-2, 3)]: atom_labels.append(("g2", f, m))
for f, m in [(0, m) for m in range(0, 1)]: atom_labels.append(("e0", f, m))
for f, m in [(1, m) for m in range(-1, 2)]: atom_labels.append(("e1", f, m))
for f, m in [(2, m) for m in range(-2, 3)]: atom_labels.append(("e2", f, m))
for f, m in [(3, m) for m in range(-3, 4)]: atom_labels.append(("e3", f, m))
N_atom = len(atom_labels)

# Target sublevels of the closed Lambda + repump scheme (Fig. 1 of the thesis)
PROBE_GROUND    = ("g2", 2, -1)   # |F=2,  m=-1>
COUPLING_GROUND = ("g1", 1, -1)   # |F=1,  m=-1>
EXCITED_TARGET  = ("e2", 2, -2)   # |F'=2, m=-2>
REPUMP_GROUND   = ("g2", 2, -2)   # |F=2,  m=-2>  (leak state recovered by the repumper)

# ==========================================
# --- HYPERFINE DIPOLE MATRIX ELEMENTS (Eq. dipole_matrix_multilevel) ---
# ==========================================
# 87Rb D2 line: J_g = 1/2 (5S1/2), J_e = 3/2 (5P3/2), I = 3/2.
# FIX: the previous implementation used bare Clebsch-Gordan coefficients
# (qutip.clebsch) for every hyperfine manifold, which gives the correct
# *relative* strengths of Zeeman sublevels WITHIN a given (F,F') pair, but
# ignores the relative strength BETWEEN different (F,F') pairs. Verified
# numerically: summing |C|^2 for all decay channels out of a given F'
# gives 1.0*gamma for F'=0,3 but 2.0*gamma for F'=1,2 -- i.e. states with
# F'=1,2 were made to decay twice as fast as F'=0,3, which is unphysical
# (all hyperfine sublevels of the D2 line share the same natural
# linewidth). The fix below implements the full formula of the thesis,
# Eq. (dipole_matrix_multilevel):
#     d_eg^(q) = <Fe||d||Fg> (-1)^(Fe-me) * 3j(Fe,1,Fg; -me,q,mg)
# with the hyperfine reduced matrix element
#     <Fe||d||Fg> = <Je||d||Jg> (-1)^(Fe+Jg+1+I) sqrt((2Fe+1)(2Jg+1)) * 6j(Jg,Je,1;Fe,Fg,I)
# The overall fine-structure reduced matrix element <Je||d||Jg]> is a
# common constant factor for every hyperfine transition of the D2 line and
# is set to 1 here: only *relative* dipole strengths are needed, both for
# the normalized branching fractions (Eq. branching_fraction) and for the
# laser coupling terms (normalized per Eq. general_laser_interaction).
_I = Rational(3, 2)   # nuclear spin, 87Rb
_Jg = Rational(1, 2)  # 5S1/2
_Je = Rational(3, 2)  # 5P3/2

_reduced_cache = {}


def _reduced_Fe_Fg(Fe, Fg):
    """<Fe||d||Fg> / <Je||d||Jg>, from the hyperfine reduced matrix element."""
    key = (Fe, Fg)
    if key in _reduced_cache:
        return _reduced_cache[key]
    Fe_r, Fg_r = Rational(Fe), Rational(Fg)
    six = wigner_6j(_Jg, _Je, 1, Fe_r, Fg_r, _I)
    if six == 0:
        val = 0.0
    else:
        prefac = (-1) ** (Fe_r + _Jg + 1 + _I) * sqrt_sym((2 * Fe_r + 1) * (2 * _Jg + 1))
        val = float(N(prefac * six))
    _reduced_cache[key] = val
    return val


_dipole_cache = {}


def dipole_element(Fg, mg, Fe, me, q):
    """
    Normalized hyperfine dipole matrix element d_eg^(q), Eq. (dipole_matrix_multilevel):
        d_eg^(q) = <Fe||d||Fg> (-1)^(Fe-me) * threej(Fe,1,Fg; -me,q,mg)
    Correctly weights different (Fg,Fe) hyperfine manifolds relative to
    each other (via the Wigner-6j factor), unlike a bare Clebsch-Gordan
    coefficient which only captures the Zeeman (m-dependent) part.
    """
    key = (Fg, mg, Fe, me, q)
    if key in _dipole_cache:
        return _dipole_cache[key]
    if mg + q != me or abs(mg) > Fg or abs(me) > Fe:
        _dipole_cache[key] = 0.0
        return 0.0
    red = _reduced_Fe_Fg(Fe, Fg)
    if red == 0.0:
        _dipole_cache[key] = 0.0
        return 0.0
    Fe_r, Fg_r = Rational(Fe), Rational(Fg)
    tj = wigner_3j(Fe_r, 1, Fg_r, -me, q, mg)
    if tj == 0:
        _dipole_cache[key] = 0.0
        return 0.0
    phase = (-1) ** (Fe - me)
    val = float(phase * red * N(tj))
    _dipole_cache[key] = val
    return val


def branching_fraction(Fg, mg, Fe, me):
    """
    Normalized branching fraction b_eg, Eq. (branching_fraction):
        b_eg = |d_eg^(q)|^2 / sum_{g',q} |d_eg'^(q)|^2 ,   sum_g b_eg = 1
    which guarantees, by construction, that
        sum_g L_{e->g}^dagger L_{e->g} = gamma_e |e><e|
    (Eq. after branching_fraction in the thesis) for L_{e->g} =
    sqrt(gamma_e * b_eg) |g><e|, regardless of any residual F'-dependence
    of the raw (unnormalized) hyperfine dipole strengths.
    """
    q = me - mg
    if abs(q) > 1:
        return 0.0
    num = dipole_element(Fg, mg, Fe, me, q) ** 2
    if num == 0.0:
        return 0.0
    denom = 0.0
    for label, Fgp, mgp in atom_labels:
        if label not in ("g1", "g2"):
            continue
        qp = me - mgp
        if abs(qp) > 1:
            continue
        denom += dipole_element(Fgp, mgp, Fe, me, qp) ** 2
    if denom == 0.0:
        return 0.0
    return num / denom


# Kept only for backward compatibility with any external script that still
# imports the old bare-Clebsch-Gordan helper. No longer used internally by
# fano.py or simulation_n_Rb_montecarlo.py (see dipole_element/branching_fraction
# above for the physically-correct, 6j-weighted replacement).
def safe_clebsch(j1, j2, j3, m1, m2, m3):
    from qutip import clebsch
    if not (abs(j1 - j2) <= j3 <= (j1 + j2)):
        return 0.0
    if m1 + m2 != m3:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    return float(clebsch(j1, j2, j3, m1, m2, m3))
