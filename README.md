# Ket Cat

**Ab Initio Neutral Atom Quantum Computer Emulator**  
|😾⟩, pronounced as "Ket Cat" is an independent quantum mechanics framework in modern C++ for simulating how neutral atom quantum circuits emerge from real atomic physics.
Instead of idealized gate algebra alone, the simulator models laser-atom interactions, Hamiltonian dynamics and wavefunction evolution directly through numerical solutions of the Time-Dependent Schrödinger Equation (TDSE).
The project combines quantum control, neutral atom physics, compile-time software architecture and scientific visualization into a single experimental framework.

[![3-qubit Quantum Fair Dice](doc/dice.gif)](https://raw.githubusercontent.com/bercesadam/KetCat/master/doc/dice.mp4)

My 3-qubit Quantum Fair Dice circuit aiming a uniform distribution among the first 6 states, showcasing the new phase disk visu and precise raman rotations of arbitrary theta angles.  (Please click on the GIF to view/download the video with full time and colour resolution.)

[![2-qubit Grover Search demonstration](doc/grover_preview.gif)](https://raw.githubusercontent.com/bercesadam/KetCat/master/doc/grover.mp4)

A successful test of 2-qubit Grover Search demonstration on Cesium atoms, performed purely by solving the Time-Dependent Schrödinger Equation, integrated on a bit more than 600 million timesteps. (Please click on the GIF to view/download the video with full time and colour resolution.)

---

### Concept: The "Digital Quantum Observatory"
The project is intended to be a bridge between Quantum Circuits and Atomic Physics, built with Software Engineering precision.
KetCat is not a mass-market research tool; it is an independent, one-man research and engineering project, meant to be a work of technological art, an architectural experiment and a tool for personal learning and explorations in quantum mechanics. While most quantum simulators stop at gate-level matrix multiplications, KetCat digs down to the "silicon" of the universe: it simulates the dynamics of **laser-atom interactions**. Here, quantum gates are not abstract unitary operators but the result of real-time physical processes (e.g., STIRAP protocols) governed by fundamental laws.

The project began as a fully `constexpr` logical circuit simulator. Today, the focus has shifted entirely to the Physical Layer:

*   **Ab Initio Simulation**: Numerically solving the TDSE with high temporal resolution (millions of time steps).
*   **Laser-Atom Interactions**: High-fidelity modeling of alkali atom (e.g., Cesium) pulse shaping and population transfer.
*   **Gate Design**: Realizing quantum gates via Raman lasers on Ladder Systems (for Bloch rotations and Rydberg excitation for CPHASE gates)
*   **Hybrid Basis Sets**: Combining **Slater-Type Orbitals (STO)** for core-electron shielding with **Quantum Defect Theory (QDT)** for high-lying Rydberg states.

### The "Architect" Approach (Engineering Principles)
As I am working as a System and SW Architect in the automotive industry, I've tried to bring my mindset into this project as well:

*   **Type Safety & Compile-Time Verification**: Utilizing C++20 Concepts and Templates to enforce eg. Hilbert space dimensions and operator compatibility at compile time.
*   **Data Visualization**: A custom, phase-encoded wave function renderer that transforms state vectors into visual aesthetics and "simple" debugging-
*   **Clean Architecture**: Separating the mathematical primitives, linear algebra, atomic physics, laser control and logical quantum circuits (and more). Te program realizes a clean pipeline which compiles logical gates into physical instructions, then laser pulses, which finally results in a Hamiltonian used for TDSE evolution. (See the architecture section below.)

---

### Quickstart

As simple as:

```bash
mkdir build && cd build
cmake ..
cmake --build .
./smoke_test
```

The project has no dependencies, other then the standard library, you only need a C++23 compliant complier.

Each simulation yields a KWF file, which is my own binary file format for simulation output data which is consumed by the supplied Python-based visu tool,
which generates a series of PNG frames (see the visu showcase on the top of this page).

---

### Architecture & Main Pipeline

The project aims to establish a clean, maintainable, and reusable architecture, where the mathematical, physical and logical layers, as well the control logic is clearly separated.

The **initialization pipeline** (executed partially at compile time, limited by `constexpr` depth of the compiler) operates as follows:

- **Atom Configuration Input**  
  The system receives an atom configuration (e.g., via the `main()` function), describing:
  - The chemical element (currently elements from Main Group 1 are supported)
  - The operation space: a set of eigenstates (e.g., 6s, 6p, 7s)
  - The mapping of these states to logical levels \(|0\rangle\), \(|1\rangle\), and \(|r\rangle\)

- **Manifold Initialization**  
  A `NeutralAtomManifold` is constructed from the configuration.  
  The full 2D eigenfunctions are calculated in an effectively infinite spatial Hilbert space, sampled on a predefined discretization grid according to the alkali atom model (see *Model Selection Logic* for details). The constructed bases set is then orthonormalized using Modified Gram-Schmidt method for correct operation space construction.

- **Projection to Reduced Hilbert Space**  
  The system projects the physical atom state into a reduced, finite-dimensional Hilbert space (operation space), where:
  - Each eigenstate is represented by a single complex amplitude
  - All TDSE (Time-Dependent Schrödinger Equation) evolution is performed

  KetCat uses a **global full state vector** (similar to Qiskit), but defined over the operation space.  
  Its size depends on the number of modeled eigenstates; however, at minimum (including logical, Rydberg, and intermediate states), it scales as: 5^QubitCount
  The logical qubit state vector is obtained as a projection from this global state.

- **Physical Parameter Computation**  
  Dipole matrix elements and Hartree energies are computed for each eigenstate and injected into the reduced model, ensuring physical realism throughout the simulation.


The **main execution flow**, governed by the `QuantumProcessor` class, is structured as follows:

- A `QuantumCircuit` (composed of `QuantumGate`s) is specified using a compile-time, type-safe DSL interface.
- Logical quantum gates are compiled into physical control instructions:
- Laser parameters and pulse sequences are generated  
- A time-dependent Hamiltonian is constructed for each time step based on the pulse envelope
- The perturbations of the state vector and Hamiltonians are being splitted into an Interaction (Dirac) picture
- The TDSE is numerically integrated over a selected subspace of the global state vector
- The "easy part", the atomic eigenenergies are solved analytically to get the Schrödinger picture for outputs calculation
- The resulting state is projected into various subspaces (ie. logical probabilities, or density matrices to get single qubit purity and basis state amplitudes for visualization and evaluation
---

### Mathematical and Physical Foundations

#### **1. Atomic Structure & Hybrid Model Selection**
For alkali atoms, the valence electron experiences a Coulomb-like potential at long range, but the tightly bound inner-shell electrons cause significant core polarization and screening. KetCat utilizes an automated, compile-time **Effective Radial Orbital** meta-generator (`EffectiveRadialOrbital` class) to select and construct the optimal **seed wavefunction**:

* **Model Selection Logic**:
    * **Hydrogenic / QDT Seed**: Selected for high angular momentum states ($l \geq 3$) or states with negligible quantum defects ($\delta_l < 0.05$). These are generated via the `HydrogenOrbitalRadial` class since the valence electron's probability density lies mostly outside the ionic core.
    * **Slater-Type Orbital (STO) Seed**: Selected for low-l states where core penetration is highly significant, managed via the `SlaterOrbitalRadial` class.

* **Effective Nuclear Charge & Screening**: For STO seeds, the framework computes the screened nuclear charge ($Z_{\text{eff}}$) at compile time using **Slater's Rules**:

    $$Z_{\text{eff}} = Z - \sigma$$

    where $Z$ is the atomic number and $\sigma$ is the shielding constant derived from the atom's specific electron configuration. The electron configuration for each alkali metal is fully resolved during compilation based solely on the $Z$ constant input, minimizing hardcoded values. The resulting radial decay parameter $\zeta$ is defined as:

    $$\zeta = \frac{Z_{\text{eff}}}{n^*}$$

* **Quantum Defects and Hydrogenic Orbitals**: For higher (eg. Rydberg) states, KetCat integrates a static lookup table calibrated against **NIST spectroscopic data** to determine the effective principal quantum number:

    $$n^* = n - \delta_l$$

    To calculate the radial part of the wavefunction for non-integer $n^*$ values, the textbook associated Laguerre polynomials are generalized using **Kummer's Confluent Hypergeometric Function** (${}_1F_1$):

```math
u_{n^*l}(r) \propto r^{l+1} \cdot e^{-\frac{r}{n^* a_{\text{eff}}}} \cdot {}_{1}F_{1}\left(-(n^* - l - 1), 2l + 2, \frac{2r}{n^* a_{\text{eff}}}\right)
```

* **2D Wavefunction Slice Generation**: The angular part of the wavefunction is universal across models and is governed by the Spherical Harmonics $Y_l^m(\theta, \phi)$, computed using Associated Legendre Polynomials $P_l^m(\cos\theta)$. The main `Hydrogenic2D` class creates planar cross-sections on an $(x, z)$ grid. To avoids numerical singularity (division by zero at the nucleus), the slice is evaluated at a tightly constrained offset $y = R_{\text{min}}$.

---

#### **2. Quantum Control & Laser-Atom Interaction**
Logical qubit states and quantum gate operations are mapped to coherent population transfers within an $N$-level ladder manifold (typically a 3-level system: $|0\rangle \leftrightarrow |1\rangle \leftrightarrow |2\rangle$).

* **STIRAP Protocol**: Coherent population transfer is driven via **Stimulated Raman Adiabatic Passage** (STIRAP). By applying a counter-intuitive pulse sequence—where the Stokes laser ($\Omega_S$) precedes the Pump laser ($\Omega_P$)—the system is trapped in a time-dependent, radiationless **dark state**:

    $$|\text{dark}(t)\rangle = \cos\vartheta(t)|0\rangle - \sin\vartheta(t)|2\rangle \quad \text{where} \quad \tan\vartheta(t) = \frac{\Omega_P(t)}{\Omega_S(t)}$$

    This approach achieves high-fidelity population transfer while bypassing the lossy intermediate state $|1\rangle$.

* **RWA Hamiltonian Construction**: The time-dependent drive is modeled in the **Rotating Wave Approximation** (RWA). The `MultiRwaRabiHamiltonian` class constructs the effective tridiagonal Hamiltonian, explicitly incorporating cumulative multi-photon detunings $\Delta_i$ and second-order **AC Stark shifts** ($\Delta E_{\text{Stark}}$) coming from off-resonant channels:

```math
\mathbf{H}_{\text{RWA}}(t) = \frac{\hbar}{2} \begin{pmatrix} 0 & \Omega_P(t) & 0 \\ \Omega_P(t) & 2\Delta_P + 2\Delta E_{\text{Stark}, 1} & \Omega_S(t) \\ 0 & \Omega_S(t) & 2(\Delta_P - \Delta_S) + 2\Delta E_{\text{Stark}, 2} \end{pmatrix}
```

---

#### **3. Numerical Propagation (Crank–Nicolson Solver)**
The real-time quantum dynamics driven by the laser pulses are simulated by integrating the Time-Dependent Schrödinger Equation (TDSE):

$$i\hbar \frac{\partial}{\partial t}|\Psi(t)\rangle = \mathbf{H}(t)|\Psi(t)\rangle$$

* **Crank–Nicolson Discretization**: To maintain strict **unitarity** (norm-preservation) without numerical drift, an implicit midpoint integration scheme is implemented:

    $$\left( \mathbf{I} + \frac{i \Delta t}{2\hbar} \mathbf{H}^{n+\frac{1}{2}} \right) |\Psi^{n+1}\rangle = \left( \mathbf{I} - \frac{i \Delta t}{2\hbar} \mathbf{H}^{n+\frac{1}{2}} \right) |\Psi^n\rangle$$

* **Solver Execution Modes**:
    * **Single-Atom Execution**: For single-qubit gates, the matrix $\mathbf{H}$ maintains a strict **tridiagonal** structure. The linear system is solved in optimal $\mathcal{O}(N)$ time using a specialized `CrankNicolsonSolver` leveraging the **Thomas Algorithm**. This routine is highly optimized and can execute millions of time steps efficiently.
    * **Multi-Atom Execution**: When tracking multi-atom (e.g., 2-qubit) configurations or multi-channel Rydberg-Rydberg interactions, the spatial block-diagonal layout introduces non-tridiagonal coupling terms. The engine dynamically switches to a pivoting **Gaussian Elimination** linear solver to preserve numerical precision across the higher-dimensional tensor-product Hilbert space.

---

### Known Limitations, Modeling Assumptions and Engineering Tradeoffs

KetCat intentionally focuses on coherent single-atom and small-system dynamics with high temporal resolution, prioritizing physical interpretability, numerical stability and architectural clarity over exhaustive physical completeness. The current implementation therefore makes several deliberate modeling assumptions and simplifications, which I consider as intentional design decisions:

* **Ladder-type coupling topology**  
  The TDSE solver exploits the tridiagonal structure of ladder-type interaction Hamiltonians for computational efficiency (see the numerical methods section). As a consequence, direct couplings are currently restricted to neighboring basis states. However, the simulator allows arbitrary user-defined basis construction, enabling the inclusion of auxiliary or weakly coupled states for leakage and population-drain modeling.

* **Assumptions on drive laser polarization**  
  Transition amplitudes are derived from dipole matrix elements computed directly from numerically integrated wavefunctions rather than from hardcoded selection rules. The present implementation performs this calculation on reduced radial wavefunctions, meaning that angular and polarization-dependent couplings are treated implicitly. Consequently, the simulator currently represents an effectively single-polarization interaction picture.

* **No hyperfine or spin-resolved structure**  
  The current quantum number abstraction and wavefunction generators model orbital states using the \((n,l,m)\) quantum numbers only. Spin, fine structure and hyperfine interactions are not yet included. While the project is conceptually inspired by neutral atom architectures such as QuEra's Cesium platforms, the present demonstrations utilize simplified orbital-state encodings (e.g. \(6s\) and \(7s\)) instead of hyperfine ground-state qubits.

* **Single-particle approximation for most operations**  
  Each qubit is primarily modeled as an individual atom evolving under its local Hamiltonian. Multi-atom effects are introduced only in dedicated interaction Hamiltonians (e.g. simplified Rydberg blockade modeling). Effects such as collective many-body dynamics, dense atomic interactions and relativistic corrections are currently outside the intended scope of the simulator.

* **No open-system decoherence or measurement collapse**  
  The framework presently models coherent unitary evolution only. Environmental decoherence, spontaneous emission channels, noise processes and wavefunction collapse during measurement are not yet incorporated. Reported state populations therefore correspond to idealized coherent probabilities derived from the propagated wavefunction.


