#pragma once
#include "quantum_bit/neutral_atom_manifold.h"
#include "laser/single_qbit_control.h"
#include "hilbert_space/basis_sets/subspace_helper.h"

namespace KetCat
{
    template <natural_t QBitCount, NeutralAtomTypeConfig Config>
    class NeutralAtomComputer
    {
        using ConfigType = std::remove_cvref_t<decltype(Config)>;

        NeutralAtomManifold<Config> m_Manifold;

		SubspaceHelper<ConfigType::LevelCount, QBitCount> m_SubSpace;

        std::array<SingleQubitControl<Config>, QBitCount> m_SingleQubitControl;

		StateVector<typename decltype(m_SubSpace)::FullHilbertSpace> m_GlobalStateVector;

    public:
        NeutralAtomComputer()
        {
            StateVector<OperationHilbertSpace> OneAtomSeed = m_Manifold.getOperationSeed();
			m_GlobalStateVector = m_SubSpace.productStateFromSeed(OneAtomSeed);
		}

        void execute(QuantumCircuit<QBitCount> circuit)
        {

        }
    };
}