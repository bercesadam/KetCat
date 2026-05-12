#pragma once
#include "quantum_bit/neutral_atom_manifold.h"
#include "laser/single_qbit_control.h"
#include "hilbert_space/basis_sets/subspace_helper.h"


namespace KetCat
{
    template <NeutralAtomTypeConfig Config>
    class NeutralAtomQubit
    {
        natrual_t QubitIndex;

        using ConfigType = std::remove_cvref_t<decltype(Config)>;

        NeutralAtomManifold<Config> m_Manifold;

        SingleQubitControl<Config> m_SingleQubitControl;


        void performInstruction(/* megkapja a gatet */)
        {
            m_SingleQubitControl.applyPulseCommand(command, Psi, (performTimeStep))
        }

        void performTimeStep(const decltype(m_SingleQubitControl)::ControlLaserArray& lasers)
        {
            static const std::array<real_t, ConfigType::LevelCount> HartreeEnergies =
                m_Manifold.getHartreeEnergies();

            static const matrix_t<ConfigType::LevelCount> DipoleMatrix =
                m_Manifold.getDipoleMatrix();

            MultiRwaRabiHamiltonian<ConfigType::LevelCount>
                Hamiltonian(
                    HartreeEnergies,
                    DipoleMatrix,
                    lasers);
            
            CrankNicolsonSolver<OperationHilbertSpace>
                    solver(
                        Hamiltonian.getMatrix(),
                        TimeMaster::Clock().getTimeStep());

			SubspaceManager::applyHamiltonian<1>(solver, { QubitIndex }, Hamiltonian);
            
        }

    };


    template <natural_t QBitCount, NeutralAtomTypeConfig Config>
    class NeutralAtomQuantumProcessor
    {
        constexpr real_t TimeStep = 50; // a.u.

        using ConfigType = std::remove_cvref_t<decltype(Config)>;
        using SubspaceManager = SubspaceHelper<ConfigType::LevelCount, QBitCount>

        std::array<NeutralAtomQubit<Config>, QBitCount> m_Qubits;

		StateVector<SubspaceManager::FullHilbertSpace> m_GlobalStateVector;

        
    public:
        NeutralAtomQuantumProcessor()
        {
            TimeMaster::Clock().init(TimeStep);
            StateVector<OperationHilbertSpace> OneAtomSeed = m_Manifold.getOperationSeed();
			m_GlobalStateVector = SubspaceManager::productStateFromSeed(OneAtomSeed);
		}

        void execute(QuantumCircuit<QBitCount> circuit)
        {

        }
    };
}