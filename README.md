# EIT Cooling of a Trapped Atom

## Status: work in progess(something may be wrong)
---

## Overview

This project implements a numerical simulation of **Electromagnetically Induced Transparency (EIT) cooling** using the QuTiP library.  
The system models a **three-level atom in $\Lambda$ configuration coupled to a quantized vibrational mode** (harmonic trap).

---

## Physical Model

We consider:
- Two ground states: $|g_1\rangle, |g_2\rangle$
- One excited state: $|e\rangle$
- A quantized harmonic oscillator (vibrational motion of the trapped atom)

The goal of EIT cooling is to **suppress carrier excitation and enhance red-sideband transitions**, allowing the system to relax toward low vibrational states.

In this regime ($\Delta > 0$), the lasers are detuned to the blue of the atomic transition.
- **Dark State:** Interference prevents absorption at the carrier frequency ($\Delta_p = \Delta_c$)
- **Fano Resonance:** A narrow absorption peak is generated

We tune the **Stark shift**:
$$\delta = \frac{\Omega_c^2}{4\Delta_c}$$

so that this peak aligns with the **red sideband** ($\omega_{\text{trap}}$), removing phonons from the system until the ground state is reached.

---

## 📂 Project Structure

`Fano_profile.py` is a simple simulation of a Fano spectrum ignoring the harmonic oscillator.

The code inside the folder `EIT_cooling` is modularized to allow fast visualization without repeating calculations:

| File | Function |
| :--- | :--- |
| `config.py` | Contains all physical parameters ($\Omega, \Delta, \eta, \nu$) and simulation time-steps |
| `simulation.py` | Builds the Hamiltonian, solves the Master Equation (Lindblad), and saves results in `results/` |
| `plot.py` | Loads saved data and generates a 3-panel animation of the cooling process |
| `requirements.txt` | Dependencies: `numpy`, `scipy`, `matplotlib`, `qutip` |
| `results` | results |
| `plots` | plots |

---

## Workflow

### Step 1 — Define Physical Parameters

We define all relevant energy scales in units of the spontaneous emission rate $\gamma$:

- $\gamma$: spontaneous emission rate  
- $\nu$: trap frequency  
- $\Delta_c$: control laser detuning  
- $\Delta_p$: probe detuning  
- $\Omega_c$: control Rabi frequency  
- $\Omega_p$: probe Rabi frequency  
- $\eta$: Lamb-Dicke parameter  

The control field is tuned to satisfy the **EIT cooling condition (Stark-shift matching)**:
$$
\Omega_c \approx \sqrt{4|\Delta_c|\nu}
$$

#### Parameters Used

| Parameter | Symbol | Value |
| :--- | :--- | :--- |
| Trap Frequency | ν | 0.5 γ |
| Coupling Detuning | Δ_c | 15.0 γ |
| Lamb-Dicke Parameter | η | 0.35 |
| Max Phonons | N_vib | 25 |

---

### Step 2 — Initial State

The system starts in a product state:

- Atom in ground state $$|g_1\rangle$$
- Motion in a Fock state $$|n=15\rangle$$

$$
|\psi_0\rangle = |g_1\rangle \otimes |n=15\rangle
$$

This allows us to observe vibrational relaxation during cooling.

---

### Step 3 — Build the Hilbert Space

We construct the full space as:

$$
\mathcal{H} = \mathcal{H}_{atom} \otimes \mathcal{H}_{motion}
$$

Using QuTiP:
- Atomic operators: projectors and transitions between levels
- Motional operators: creation and annihilation operators

Key operators:
- $a, a^\dagger$: phonon ladder operators  
- $n = a^\dagger a$: phonon number operator  
- $\sigma_{e g_i}$: atomic transitions  

---

### Step 4 — Hamiltonian

The full Hamiltonian includes:

#### 1. Free energy terms
$$
H_0 = \Delta_p |g_1\rangle\langle g_1| + \Delta_c |g_2\rangle\langle g_2| + \nu a^\dagger a
$$

#### 2. Laser interactions

- Control field couples $|g_2\rangle \leftrightarrow |e\rangle$
- Probe field couples $|g_1\rangle \leftrightarrow |e\rangle$ including motional coupling:

$$
H_{\text{int}} \propto \Omega_p \, \sigma_{e g_1}
\left(1 + i \eta (a + a^\dagger)\right) + \text{h.c.}
$$

This term introduces **sideband transitions**, which are essential for cooling.

---

### Step 5 — Dissipation

Spontaneous emission is included via collapse operators:

- $|e\rangle \rightarrow |g_1\rangle$
- $|e\rangle \rightarrow |g_2\rangle$

This introduces irreversibility and allows cooling.

---

### Step 6 — Solve Master Equation

We solve the Lindblad master equation:

$$
\frac{d\rho}{dt} = -i[H,\rho] + \sum_k \mathcal{D}[C_k]\rho
$$

Using:
```python
mesolve(H, psi0, t_list, c_ops, [n_op])
```

### Step 7 — EIT Spectrum

To characterize the cooling mechanism, we compute the **steady-state excitation spectrum** by scanning the probe detuning $\Delta_p$.

For each value of $\Delta_p$:

1. We construct the atomic Hamiltonian in the rotating frame  
2. We solve the steady-state master equation:
$$\frac{d\rho}{dt} = 0$$
3. We evaluate the excited-state population:
$$I(\Delta_p) = \langle e|\rho_{ss}|e\rangle$$

This spectrum reveals:
- The EIT dark state (transparency window)
- Suppression of carrier excitation
- Asymmetric sideband structure (Fano-like profile)

---

### Step 8 — Physical Observables

#### Vibrational dynamics

The expectation value of the phonon number is:
$$\langle n(t) \rangle$$

This is the main indicator of cooling efficiency, showing how the system approaches a low-energy steady state.

---

#### Quantum state evolution

The full density matrix $\rho(t)$ is stored at each timestep, allowing reconstruction of:
- Atomic populations  
- Vibrational populations  
- Atom–motion correlations  

---

#### EIT spectrum

The steady-state scan provides:
- Excitation probability vs detuning  
- Identification of dark-state suppression  
- Cooling-resonant red sideband enhancement  


#### Note on Sideband Detuning in the Rotating Frame

In this simulation, we follow the convention $\Delta = \omega_L - \omega_0$. Since the system operates in the **Blue-Detuned regime** ($\Delta > 0$), the mapping of sidebands in the rotating frame follows the specific physics of the Morigi EIT scheme:

*   **Red Sideband (Cooling):** $\Delta_p = \Delta_c + \nu$  
    To extract energy ($n \to n-1$), the probe must be tuned to the "bright" dressed state. Because the AC Stark shift $\delta$ pushes this peak toward higher frequencies in a blue-detuned setup, the cooling resonance occurs at **$+\nu$** relative to the EIT dark state.

*   **Blue Sideband (Heating):** $\Delta_p = \Delta_c - \nu$  
    The transition that adds energy ($n \to n+1$) is intentionally mapped to a lower frequency in the rotating frame, falling into the EIT "dark" region where absorption is suppressed.

This "sign flip" (where $+\nu$ cools and $-\nu$ heats) is a direct consequence of the asymmetric Fano profile generated by the blue-detuned coupling laser, ensuring that the maximum scattering rate is reserved for phonon-subtraction processes.

---
## How to Run

### 1. Setup Environment

Install the required libraries:

```bash
pip install -r requirements.txt
```

### Execute Simulation
Run the numerical solver. This will generate data files in a new results/ folder:
```bash
python simulation.py
```

### Visualize Results
Run the plotting script to open the interactive animation:
```bash
python plotting.py
```

---

## Visualization Details

The animation generated by `plotting.py` in plots consists of:

- **Phonon Distribution (Left):**  
  A bar chart showing the population of Fock states $|n\rangle$ evolving from $n = 15$ towards $n = 0$.

- **EIT Spectrum (Center):**  
  The Fano profile showing the "dark" window at the carrier and the cooling peak at the Red Sideband.

- **Cooling Curve (Right):**  
  The expectation value of the phonon number $\langle n \rangle$ vs time.

---
<p align="center">
  <img src="Fano/fano_profile.png" alt="Fano profile" width="450"/>
</p>

![Cooling Evolution](EIT_cooling/plots/cooling_evolution.gif)

## ⏱️ Time Units in the Simulation

Time in this simulation is expressed in natural units set by the decay rate `gamma`. All frequencies in the Hamiltonian (such as `gamma`, `Delta_c`, `nu`, and `Omega_c`) are defined in the same units, so time is measured as:

t_unit = 1 / gamma

For a concrete physical interpretation, if `gamma = 6 MHz` (i.e. `6 × 10^6 s⁻¹`), then one unit of simulation time corresponds to:

1 time unit = 1 / gamma = 1.67 × 10⁻⁷ s ≈ 167 ns

With the parameter:

t_stop = 3500

the total simulated physical time becomes:

t_max = 3500 / gamma ≈ 0.58 ms

A few useful conversions is:
10000  → 1.67 ms
In general, the conversion between simulation time and physical time is:

t_physical = t_simulation / gamma

# 87Rb EIT Cooling Simulation: 24-Level Hyperfine Manifold

## Status: work in progress (model under refinement)

---

## Overview

This module extends the standard three-level EIT cooling model to a fully realistic description of $^{87}$Rb on the $D_2$ transition. Instead of an ideal $\Lambda$ system, the full hyperfine and Zeeman structure is included, resulting in a 24-level atomic model.

The aim is to investigate how real atomic complexity affects:
- EIT cooling efficiency  
- Fano interference  
- Optical pumping dynamics  
- Leakage into unwanted states  

---

## 1. From a 3-Level Model to a 24-Level System

### 1.1 Ideal $\Lambda$ System

The standard EIT cooling scheme is based on a three-level system:
- Two ground states: $$|g_1\rangle$, $|g_2\rangle$$
- One excited state: $$|e\rangle$$

Quantum interference creates a dark state:
$$
|D\rangle \propto \Omega_c |g_1\rangle - \Omega_p |g_2\rangle
$$

which suppresses carrier excitation and enables cooling.

---

### 1.2 Real Atomic Structure

In $$^{87}$$Rb, this picture is only approximate. The atomic structure includes:

- Ground states:
  $$F=1 \ (3 \ states), \quad F=2 \ (5 \ states)$$

- Excited states:
  $$F'=0,1,2,3 \ (\text{16 total Zeeman states})$$

Total Hilbert space:
$$\mathcal{H}_{atom} = 24 \text{ levels}$$

The effective $\Lambda$ system is embedded inside this larger structure and continuously perturbed by additional couplings.

---

## 2. Model Construction

### 2.1 Hamiltonian

The Hamiltonian is constructed in the rotating frame of the coupling laser:

$$
H = H_{\text{atom}} + H_{\text{Zeeman}} + H_{\text{int}}
$$

#### Atomic energies

Hyperfine splittings are defined relative to $F'=2$.

#### Zeeman shifts

$$E_Z = g_F \mu_B B m_F$$

#### Laser interaction

Each transition is weighted by Clebsch–Gordan coefficients:

$$
H_{\text{int}} = \frac{\Omega}{2} \sum_{g,e} C_{g,e} |e\rangle\langle g| + \text{h.c.}
$$

Laser configuration:
- Coupling field: $\sigma^-$ polarization  
- Probe field: $\pi$ polarization  
- Repumper: $\sigma^-$ polarization  

---

### 2.2 Dissipation

Spontaneous emission is described using Lindblad operators:

$$
\mathcal{L}(\rho) = \sum_{e,g,q} \gamma C_{g,e}^{(q)} \left( |g\rangle\langle e| \rho |e\rangle\langle g| - \frac{1}{2}\{ |e\rangle\langle e|, \rho \} \right)
$$

All allowed transitions are included:
- $F' \rightarrow F=1,2$
- Polarizations $q = -1, 0, +1$

---

### 2.3 Steady-State Solution

The system is solved using the Lindblad master equation:

$$
\dot{\rho} = -i[H,\rho] + \mathcal{L}(\rho)
$$

At steady state:
$$\dot{\rho} = 0$$

Numerically:
```python
rho_ss = steadystate(H, c_ops)
```

The absorption spectrum is defined as:
$$A(\Delta_p) = \text{Tr}(\rho_{ss} P_e)$$

---

## 3. Zeeman Structure and Magnetic Field

### 3.1 Importance of Zeeman Sublevels

Zeeman states introduce:
- Multiple coupling pathways  
- Additional dark states  
- Optical pumping channels  

The ideal $\Lambda$ system is no longer isolated.

---

### 3.2 Two-Photon Resonance

The resonance condition becomes:
$$\Delta_p - \Delta_c + \Delta_Z = 0$$

where $\Delta_Z$ is the Zeeman shift difference.

---

### 3.3 Case $B = 0$

If $B=0$:
- All Zeeman states are degenerate  
- Multiple dark states appear  
- Population spreads across sublevels  

Consequences:
- Reduced EIT contrast  
- Increased leakage  
- Lower cooling efficiency  

---

## 4. Role of the Repumper

Spontaneous emission transfers atoms to $F=1$:
$$|e\rangle \rightarrow |F=1\rangle$$

Without repumping:
- Population accumulates in $F=1$  
- Cooling stops  

The repumper must be properly polarized:

- Correct polarization → efficient recycling  
- Incorrect polarization → dark states  

This determines whether the cooling cycle remains closed.

---

## 5. Trap Frequency and Motional Effects

### 5.1 Sideband Structure

The trap frequency $\nu$ defines:
- Red sideband: $\Delta = +\nu$  
- Blue sideband: $\Delta = -\nu$  

EIT cooling condition:
$$\delta = \frac{\Omega_c^2}{4\Delta_c} \approx \nu$$

---

### 5.2 Leakage Mechanisms

For larger $\nu$:
- Stronger coupling to sidebands  
- Increased spectral overlap  

This enhances excitation of unwanted states, especially $F'=3$, since the Clebsch–Gordan Scaling is

Effective coupling:
$$\Omega_{\text{eff}} = \Omega \cdot C_{g,e}$$

WRITE BETTER

---

## 6. Repository Structure

The folder `EIT_cooling_Rb` contains:

| File | Description |
|------|------------|
| `config.py` | Defines atomic levels, parameters, Clebsch coefficients |
| `fano.py` | Builds Hamiltonian, solves steady state, saves data |
| `plot_fano.py` | Loads data and generates plots |
| `level_diagram.py` | Generates level diagram image |
| `results_fano/` | Numerical outputs |
| `images/` | Generated plots |

---

## 7. How to Run

### Install dependencies
```bash
pip install -r requirements.txt
```

### Run simulation
Update the name of the experiment in config.py then run

```bash
python fano.py
```

### Plot results
```bash
python plot_fano.py
```

### Generate level diagram
```bash
python level_diagram.py
```

---

## 8. Results and Visualization

### 8.1 Absorption Spectrum

$$A(\Delta_p) = \langle P_e \rangle$$

Features:
- EIT transparency window  
- Fano asymmetric profile  
- Motional sidebands  

---

### 8.2 Ground State Populations

Tracks:
- Optical pumping efficiency  
- Population in unwanted $m_F$ states  
- Repumper effectiveness  

---

### 8.3 Leakage Ratio

Defined as:
$$\text{Leakage} = \frac{P_{F'=3}}{P_{\text{total}}}$$

Measures excitation lost to off-resonant states.


<p align="center">
  <img src="EIT_cooling_Rb/images/eit_diagram.png" width="900"/>
</p>

<p align="center">
  <img src="EIT_cooling_Rb/images/fano_images/plot_fano_profile.png" width="900"/>
</p>

