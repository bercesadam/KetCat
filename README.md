# Ket Cat

**Ab Initio Neutral Atom Quantum Computer Emulator**  
|😾⟩, pronounced as "Ket Cat" is an independent quantum mechanics framework in modern C++ for simulating how neutral atom quantum circuits emerge from real atomic physics.
Instead of idealized gate algebra alone, the simulator models laser-atom interactions, Hamiltonian dynamics and wavefunction evolution directly through numerical solutions of the Time-Dependent Schrödinger Equation (TDSE).
The project combines quantum control, neutral atom physics, compile-time software architecture and scientific visualization into a single experimental framework.

<img src="https://raw.githubusercontent.com/bercesadam/QuantumCircuitsinCompiler/master/doc/demo.gif" alt="One qubit demonstration" width="1024" style="text-align:center">
A successful test of single-qubit gates on a Cesium atom with STIRAP Laser drive, performed purely with solving the Time-Dependent Schrödinger (actually integrated on ~80 million time steps).

---

### Concept: The "Digital Quantum Observatory"
The project is intended to be a bridge between Quantum Circuits and Atomic Physics, built with Software Engineering precision.
KetCat is not a mass-market research tool; it is an independent, one-man research and engineering project, meant to be a work of technological art, an architectural experiment and a tool for personal learning and explorations in quantum mechanics. While most quantum simulators stop at gate-level matrix multiplications, KetCat digs down to the "silicon" of the universe: it simulates the dynamics of **laser-atom interactions**. Here, quantum gates are not abstract unitary operators but the result of real-time physical processes (e.g., STIRAP protocols) governed by fundamental laws.


### Current Focus: Neutral Atoms & STIRAP
The project began as a fully `constexpr` logical circuit simulator. Today, the focus has shifted entirely to the Physical Layer:

*   **Ab Initio Simulation**: Numerically solving the TDSE with high temporal resolution (millions of time steps).
*   **Laser-Atom Interactions**: High-fidelity modeling of alkali atom (e.g., Cesium) pulse shaping and population transfer.
*   **Gate Design**: Realizing quantum gates (Hadamard, Pauli rotations) via fractional and inverted **STIRAP** (Stimulated Raman Adiabatic Passage) protocols.
*   **Hybrid Basis Sets**: Combining **Slater-Type Orbitals (STO)** for core-electron shielding with **Quantum Defect Theory (QDT)** for high-lying Rydberg states.

### The "Architect" Approach (Engineering Principles)
As I am working as a System and SW Architect in the automotive industry, I've tried to bring my mindset into this project as well:

*   **Type Safety & Compile-Time Verification**: Utilizing C++20 Concepts and Templates to enforce eg. Hilbert space dimensions and operator compatibility at compile time.
*   **Data Visualization**: A custom, phase-encoded wave function renderer that transforms abstract complex numbers into visual aesthetics and "simple" debugging.
*   **Clean Architecture**: Separating the mathematical primitives, linear algebra, atomic physics, laser control and logical quantum circuits (and more).

### Known Limitations, Modeling Assumptions and Engineering Tradeoffs

KetCat intentionally focuses on coherent single-atom and small-system dynamics with high temporal resolution, prioritizing physical interpretability, numerical stability and architectural clarity over exhaustive physical completeness. The current implementation therefore makes several deliberate modeling assumptions and simplifications:

* **Ladder-type coupling topology**  
  The TDSE solver exploits the tridiagonal structure of ladder-type interaction Hamiltonians for computational efficiency (see the numerical methods section). As a consequence, direct couplings are currently restricted to neighboring basis states. However, the simulator allows arbitrary user-defined basis construction, enabling the inclusion of auxiliary or weakly coupled states for leakage and population-drain modeling.

* **Polarization simplification**  
  Transition amplitudes are derived from dipole matrix elements computed directly from numerically integrated wavefunctions rather than from hardcoded selection rules. The present implementation performs this calculation on reduced radial wavefunctions, meaning that angular and polarization-dependent couplings are treated implicitly. Consequently, the simulator currently represents an effectively single-polarization interaction picture.

* **No hyperfine or spin-resolved structure**  
  The current quantum number abstraction and wavefunction generators model orbital states using the \((n,l,m)\) quantum numbers only. Spin, fine structure and hyperfine interactions are not yet included. While the project is conceptually inspired by neutral atom architectures such as QuEra's Cesium platforms, the present demonstrations utilize simplified orbital-state encodings (e.g. \(6s\) and \(7s\)) instead of hyperfine ground-state qubits.

* **Single-particle approximation for most operations**  
  Each qubit is primarily modeled as an individual atom evolving under its local Hamiltonian. Multi-atom effects are introduced only in dedicated interaction Hamiltonians (e.g. simplified Rydberg blockade modeling). Effects such as collective many-body dynamics, dense atomic interactions and relativistic corrections are currently outside the intended scope of the simulator.

* **No open-system decoherence or measurement collapse**  
  The framework presently models coherent unitary evolution only. Environmental decoherence, spontaneous emission channels, noise processes and wavefunction collapse during measurement are not yet incorporated. Reported state populations therefore correspond to idealized coherent probabilities derived from the propagated wavefunction.
---

### Quickstart

As simple as:

```bash
mkdir build && cd build
cmake ..
cmake --build .
./smoke_test
```

The project has no dependencies, other the standard library, you only need a C++23 compliant complier.

Each simulation yields a KWF file, which is my own binary file format for simulation output data which is consumed by the supplied Python-based visu tool,
which generates a series of PNG frames (see the visu showcase on the top of this page).

---

### Mathematical and Physical Foundations

KetCat operates across multiple layers of abstraction to simulate a quantum processor from first principles.

#### **1. Atomic Structure & Hybrid Model Selection**
For alkali atoms, the valence electron experiences a Coulomb-like potential at long range, but inner electrons cause significant screening. KetCat utilizes an **Effective Radial Orbital** meta-generator to create the initial **seed wavefunction**:

*   **Model Selection Logic**:
    *   **Hydrogenic/QDT Seed**: Selected for high angular momentum states ($l \geq 3$) or states with negligible quantum defects ($\delta_l < 0.05$).
    *   **Slater-Type Orbital (STO) Seed**: Selected for low-$l$ states where core penetration is significant.
*   **Effective Nuclear Charge**: For STO seeds, the simulator calculates the screened nuclear charge ($Z_{eff}$) using **Slater's Rules**:
    $$Z_{eff} = Z - \sigma$$
    where $\sigma$ is the shielding constant derived from the electron configuration. The resulting orbital decay is governed by $\zeta = Z_{eff} / n^*$.

#### **2. Wavefunction Physics (Seed Wavefunction)**
The spatial wavefunction is decomposed into radial and angular components:
$$\Psi_{n^*lm}(r, \theta, \phi) = \frac{u_{n^*l}(r)}{r} Y_l^m(\theta, \phi)$$
*   **Angular Part**: Computed via Associated Legendre polynomials $P_l^m(x)$ normalized into Spherical Harmonics $Y_l^m$.
*   **Radial Part ($u$)**: To support non-integer principal quantum numbers ($n^* = n - \delta_l$), KetCat generalizes the radial solution using the **Kummer Confluent Hypergeometric function** ($_{1}F_{1}$):
    $$u_{n^*l}(r) \propto r^{l+1} \cdot e^{-\frac{r}{n^* a_{eff}}} \cdot {}_{1}F_{1}(-(n^* - l - 1), 2l + 2, \frac{2r}{n^* a_{eff}})$$
    This ensures that nodal structures and phases correctly reflect core penetration effects.

#### **3. Quantum Control & Laser Interaction**
Quantum gates are implemented as physical transitions in a 3-level system ($|0\rangle \leftrightarrow |1\rangle \leftrightarrow |2\rangle$).
*   **STIRAP Protocol**: Uses a counter-intuitive pulse sequence (Stokes before Pump) to transfer population via a "dark state," avoiding the lossy intermediate state.
*   **Hamiltonian Construction**: The interaction Hamiltonian in the Rotating Wave Approximation (RWA) is:
    $$\mathbf{H}_{int} = \frac{\hbar}{2} \begin{pmatrix} 0 & \Omega_P(t) & 0 \\ \Omega_P(t) & 2\Delta_P & \Omega_S(t) \\ 0 & \Omega_S(t) & 2(\Delta_P - \Delta_S) \end{pmatrix}$$

#### **4. Numerical Propagation (Crank–Nicolson)**
Temporal evolution is handled by solving the TDSE using the **Crank–Nicolson method**, ensuring **unitarity** (norm-preservation):
$$\left( \mathbf{I} + \frac{i \Delta t}{2\hbar} \mathbf{H} \right) \Psi^{n+1} = \left( \mathbf{I} - \frac{i \Delta t}{2\hbar} \mathbf{H} \right) \Psi^n$$
By utilizing a **tridiagonal Hamiltonian** matrix, the system is solved in $\mathcal{O}(N)$ time using the **Thomas algorithm**, allowing for millions of high-resolution time steps.



