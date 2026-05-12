# Ket Cat

**Ab Initio Neutral Atom Quantum Computer Emulator**  
|😾⟩, pronounced as "Ket Cat" is a first-principles simulation framework for modeling the coherent dynamics of neutral atom qubits via TDSE, laser-atom interaction, and adiabatic passage protocols.

<img src="https://raw.githubusercontent.com/bercesadam/QuantumCircuitsinCompiler/master/doc/demo.gif" alt="One qubit demonstration" width="1024" style="text-align:center">
A successful test of a single-qubit gates on a Cesium atom with STIRAP Laser drive, performed purely with solving the Time-Dependent Schrödinger (actually integrated on ~80 million time steps).

---

### A bridge between Quantum Circuits and Atomic Physics, built with Software Engineering precision.

KetCat is a modern C++ framework designed to unify the logical abstractions of quantum computing with the underlying physical reality of the **Time-Dependent Schrödinger Equation (TDSE)**. The project focuses on the ab initio modeling and visualization of neutral atom quantum processors. 

### 🌌 Concept: The "Digital Quantum Observatory"
KetCat is not a mass-market research tool; it is a work of technological art and an architectural experiment. While most quantum simulators stop at gate-level matrix multiplications, KetCat digs down to the "silicon" of the universe: it simulates the dynamics of **laser-atom interactions**. Here, quantum gates are not abstract unitary operators but the result of real-time physical processes (e.g., STIRAP protocols) governed by fundamental laws.


### 🚀 Current Focus: Neutral Atoms & STIRAP
The project began as a fully `constexpr` logical circuit simulator. Today, the focus has shifted entirely to the Physical Layer:

*   **Ab Initio Simulation**: Numerically solving the TDSE with high temporal resolution (millions of time steps).
*   **Laser-Atom Interactions**: High-fidelity modeling of alkali atom (e.g., Cesium) pulse shaping and population transfer.
*   **Gate Design**: Realizing quantum gates (Hadamard, Pauli rotations) via fractional and inverted **STIRAP** (Stimulated Raman Adiabatic Passage) protocols.
*   **Hybrid Basis Sets**: Combining **Slater-Type Orbitals (STO)** for core-electron shielding with **Quantum Defect Theory (QDT)** for high-lying Rydberg states.

### 🛠 The "Architect" Approach (Engineering Principles)
As I am working as a System and SW Architect in the automotive industry, I've tried to bring my mindset into this project as well:

*   **Type Safety & Compile-Time Verification**: Utilizing C++20 Concepts and Templates to enforce eg. Hilbert space dimensions and operator compatibility at compile time.
*   **Data Visualization**: A custom, phase-encoded wave function renderer that transforms abstract complex numbers into visual aesthetics and "simple" debugging.
*   **Clean Architecture**: Separating the mathematical primitives, linear algebra, atomic physics, laser control and logical quantum circuits (and more).
---

### 📦 Quickstart
```bash
mkdir build && cd build
cmake ..
cmake --build .

---

### 📚 Mathematical and Physical Foundations

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



