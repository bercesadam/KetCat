#pragma once
#include "quantum_bit/neutral_atom_manifold.h"
#include "laser/single_qbit_control.h"
#include "hilbert_space/basis_sets/subspace_helper.h"


namespace KetCat
{
    template <NeutralAtomTypeConfig Config>
    class NeutralAtomQubit
    {
        using ConfigType = std::remove_cvref_t<decltype(Config)>;

        using SubspaceManager = SubspaceHelper<ConfigType::LevelCount, QBitCount>

        NeutralAtomManifold<Config> m_Manifold;

        SingleQubitControl<Config> m_SingleQubitControl;


        void performInstruction()
        {
           
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

            psi = solver(psi);
            
        }

    };


    template <natural_t QBitCount, NeutralAtomTypeConfig Config>
    class NeutralAtomComputer
    {
        constexpr real_t TimeStep = 50; // a.u.

        using ConfigType = std::remove_cvref_t<decltype(Config)>;
        using SubspaceManager = SubspaceHelper<ConfigType::LevelCount, QBitCount>

        NeutralAtomManifold<Config> m_Manifold;

        std::array<SingleQubitControl<Config>, QBitCount> m_SingleQubitControl;

		StateVector<SubspaceManager::FullHilbertSpace> m_GlobalStateVector;

        
    public:
        NeutralAtomComputer()
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