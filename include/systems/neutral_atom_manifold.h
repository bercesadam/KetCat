#pragma once
#include <memory>
#include <type_traits>

#include "neutral_atom_config.h"

#include "hilbert_space/hilbert.h"
#include "hilbert_space/gram_schmidt_orthonorm.h"
#include "hilbert_space/reduced_energy_space.h"

#include "atomic_physics_core/elements.h"
#include "atomic_physics_core/quantum_number.h"

#include "wavefunction/atomic/2d_hydrogenic.h"
#include "wavefunction/atomic/hartree.h"

#include "hamiltonian/dipole_operator.h"

namespace KetCat
{
    /// @warning EXTREME TEMPLATE ZONE — PROCEED WITH CAUTION
    ///
    /// This class is the most template-heavy region of the entire project.
    /// Its design is intentional and optimized for strong compile-time
    /// guarantees and zero runtime overhead.
    ///
    /// Contents include (but are not limited to):
    ///   - template-template parameters
    ///   - dependent names
    ///   - dependent template names
    ///   - concepts layered on top of all of the above
    ///   - non-type template parameters (NTTPs) that pretend to be types
    ///
    /// Side effects of reading or modifying this code may include:
    ///   - confusion
    ///   - sudden questions about why something does not deduce
    ///   - appreciation for runtime polymorphism
    ///   - questioning your life choices
    ///
    /// This file prioritizes:
    ///   - compile-time safety over immediate readability
    ///   - explicitness over convenience
    ///   - correctness over comfort
    ///
    /// If this code compiles, it is doing exactly what it was designed to do.
    /// Modify only if you are comfortable with advanced C++20 template mechanics.
    ///
    /// @brief Represents the operation space of a single qubit (neutral atom)
    ///        holding the eigen states of the operation space and a provider
    ///        for single atom-related properties (dipole matrix for Hamiltonians;
    ///        full 2D state vector for visualization purposes).
    ///
    /// @tparam Config
    ///   Compile-time neutral atom configuration describing element,
    ///   discretization, logical levels, and participating quantum numbers.
    ///
    template <NeutralAtomTypeConfig Config>
    class NeutralAtomManifold
    {
    private:
        /// @brief Type-level access to the configuration NTTP.
        ///
        /// @details
        /// Since `Config` is a non-type template parameter, this alias is required
        /// to access nested types, templates, and static members in a type context.
        using ConfigType = std::remove_cvref_t<decltype(Config)>;

    public:
        using SingleAtomRadialHilbertSpace = InfiniteHilbertSpace<1_D, 8192, 2000.0>;
        //typename ConfigType::template HilbertSpaceStub<1_D>;

        using SingleAtomFullHilbertSpace = InfiniteHilbertSpace<2_D, 256, 1.0>;

    private:
        /// @brief Reduced-energy Hilbert space used for time evolution.
        ///
        /// @details
        /// This space collapses the full spatial wavefunctions into a smaller
        /// Hilbert space spanned only by energy eigenstates.
        using ReducedEnergySpaceType =
            ReducedEnergySpace<
            SingleAtomRadialHilbertSpace,
            ConfigType::LevelCount
            >;

    public:
        /// @brief Hilbert space used for single-atom operations (qubit-level dynamics).
        ///
        /// @details
        /// This is the reduced Hilbert space on which control Hamiltonians
        /// and time evolution are applied.
        using SingleAtomOperationHilbertSpace =
            typename ReducedEnergySpaceType::ReducedHilbertSpace;

    private:
        //typename ConfigType::template HilbertSpaceStub<2_D>;

        /// @brief Basis-set type alias for a given spatial dimensionality.
        ///
        /// @tparam SpatialDimensions
        ///   Spatial dimensionality of the underlying Hilbert space (e.g. 1_D, 2_D).
        ///
        /// @details
        /// A basis set is represented as a fixed-size array of wavefunctions,
        /// one for each eigenstate defined in the configuration.
        //template <DimensionTag SpatialDimensions>
		template <hilbert_space_t HilbertSpace>
        using BasisSet =
            basis_set_t<
                HilbertSpace,
                ConfigType::LevelCount
            >;

        /// @brief Owning handle for a dynamically constructed basis set.
        ///
        /// @tparam SpatialDimensions
        ///   Spatial dimensionality of the basis set.
        //template <DimensionTag SpatialDimensions>
        template <hilbert_space_t HilbertSpace>
        using BasesSetHandle = std::unique_ptr<BasisSet<HilbertSpace>>;

        /// @brief Owning handle for the reduced energy space.
        using ReducedEnergySpaceHandle =
            std::unique_ptr<ReducedEnergySpaceType>;

    private:
        /// @brief Electric dipole transition matrix between eigenstates.
        ///
        /// @details
        /// This matrix is computed from the radial parts of the wavefunctions
        /// and encodes allowed optical transitions under the electric-dipole
        /// approximation.
        inline static matrix_t<ConfigType::LevelCount> m_dipoleMatrix;

		/// @brief Eigenvalues of the energy levels in Hartree atomic units.
        inline static std::array<real_t, ConfigType::LevelCount> m_hartreeEnergies;

        /// @brief Orthonormalized full spatial basis states (1D).
        ///
        /// @details
        /// Used primarily for visualization and spatially resolved analysis.
        BasesSetHandle<SingleAtomRadialHilbertSpace> m_basisStates1D;

        /// @brief Full spatial basis states (2D) for visualization.
        ///
        /// @details
        /// Used primarily for visualization and spatially resolved analysis.
        BasesSetHandle<SingleAtomFullHilbertSpace> m_basisStates2D;

        /// @brief Reduced-energy space used for the actual qubit operations.
        ReducedEnergySpaceType m_operationSpace;


        /// @brief Construct a single basis state using a specified generator.
        ///
        /// @tparam HilbertSpace
        ///   Spatial Hilbert space on which the wavefunction is defined.
        /// @tparam WaveFunctionGenerator
        ///   Wavefunction generator class template.
        /// @tparam QuantumNumber
        ///   Quantum number type n, l, [m]
        ///
        /// @return
        ///   A fully constructed wavefunction corresponding to the given quantum numbers.
        template<
            spatial_hilbert_space_t HilbertSpace,
            template<
                spatial_hilbert_space_t,
                Element
            > class WaveFunctionGenerator,
            quantum_number_t QuantumNumber
        >
        constexpr auto makeBasisState() noexcept
        {
            return WaveFunctionGenerator<
                HilbertSpace,
                ConfigType::ChemicalElement
            >{}(QuantumNumber{});
        }

        /// @brief Construct a complete basis set using a specified generator.
        ///
        /// @tparam SpatialDimensions
        ///   Spatial dimensionality of the basis set.
        /// @tparam WaveFunctionGenerator
        ///   Generator used to construct individual basis states.
        /// @tparam IndexSequence
        ///   Compile-time index sequence used to expand over quantum numbers.
        template<
            /*DimensionTag SpatialDimensions,*/
			spatial_hilbert_space_t HilbertSpace,
            template<
                spatial_hilbert_space_t,
                Element
            > class WaveFunctionGenerator,
            std::size_t... IndexSequence
        >
        constexpr auto makeBasisSet(std::index_sequence<IndexSequence...>) noexcept
        {
            return BasisSet<HilbertSpace>{{
                this->template makeBasisState<
                    HilbertSpace,
                    WaveFunctionGenerator,
                    std::tuple_element_t<
                        IndexSequence,
                        typename ConfigType::QuantumNumbers
                    >
                >()...
            }};
        }

        /// @brief Compute the radial electric dipole transition matrix.
        ///
        /// @details
        /// Uses effective radial orbitals to compute dipole matrix elements
        /// under the electric-dipole approximation.
        void buildDipleMatrix() noexcept
        {
            m_basisStates1D =
                std::make_unique<BasisSet<SingleAtomRadialHilbertSpace>>(
                    makeBasisSet<SingleAtomRadialHilbertSpace, EffectiveRadialOrbital>(
                        std::make_index_sequence<ConfigType::LevelCount>{}
                    )
                );

            m_dipoleMatrix = buildRadialDipoleMatrix(*m_basisStates1D);

            auto MGS =
                std::make_unique<Orthonormalizer<ConfigType::LevelCount>>();

            MGS->learn(*m_basisStates1D);

            m_basisStates1D = std::make_unique<BasisSet<SingleAtomRadialHilbertSpace>>(
                MGS->apply(*m_basisStates1D)
            );
        }

        /// @brief Construct and orthonormalize the full spatial basis set (2D).
        ///
        /// @details
        /// The basis is generated using hydrogenic wavefunctions and
        /// orthonormalized using a two-pass Modified Gram-Schmidt.
        void buildFullBasisSet() noexcept
        {
            m_basisStates2D =
                std::make_unique<BasisSet<SingleAtomFullHilbertSpace>>(
                    makeBasisSet<SingleAtomFullHilbertSpace, Hydrogenic2D>(
                        std::make_index_sequence<ConfigType::LevelCount>{}
                    )
                );

            auto MGS =
                std::make_unique<Orthonormalizer<ConfigType::LevelCount>>();

            MGS->learn(*m_basisStates2D);

            m_basisStates2D = std::make_unique<BasisSet<SingleAtomFullHilbertSpace>>(
                MGS->apply(*m_basisStates2D)
            );
        }

        /// @brief Construct the reduced-energy operational Hilbert space.
        ///
        /// @details
        /// This space is used for actual time evolution and control dynamics.
        void buildOperationSpace() noexcept
        {
            //m_operationSpace =
             //   std::make_unique<ReducedEnergySpaceType>(*m_basisStates);
        }
        

        /// @brief Compile-time helper to compute Hartree energies for all eigenstates.
        ///
        /// @details
        /// This function expands the list of quantum-number *types* provided by the
        /// NeutralAtomTypeConfig using a compile-time index sequence. For each quantum
        /// number, it evaluates the corresponding Hartree energy via
        /// `calculateHartreeEnergy`.
        ///
        /// @tparam IndexSequence
        ///   Compile-time indices used to expand the quantum-number tuple.
        /// @param std::index_sequence<Is...>
        ///   Tag used to trigger parameter-pack expansion.
        /// @return
        ///   Fixed-size array of Hartree energies, one per eigenstates, in Hartree
        ///   atomic units.
        template<std::size_t... IndexSequence>
        static constexpr std::array<real_t, ConfigType::LevelCount>
            calculateHartreeEnergiesImpl(std::index_sequence<IndexSequence...>) noexcept
        {
            return {
                calculateHartreeEnergy(
                    ConfigType::ChemicalElement,
                    std::tuple_element_t<
                        IndexSequence,
                        typename ConfigType::QuantumNumbers
                    >{}
                )...
            };
        }

		/// @brief Calculate Hartree energies for all eigenstates
        void calculateHartreeEnergies() noexcept
        {
            m_hartreeEnergies = calculateHartreeEnergiesImpl(
                std::make_index_sequence<ConfigType::LevelCount>{}
            );
        }

    public:
        /// @brief Retrieve the |0⟩ operation seed state.
        ///
        /// @return
        ///   State vector representing the logical |0⟩ state in the reduced space for one qubit
        StateVector<SingleAtomOperationHilbertSpace> getOperationSeed() noexcept
        {
            return m_operationSpace.project(*m_basisStates2D,
                (*m_basisStates2D)[ConfigType::Logical0Level].m_Psi);
        }

		/// @brief Retrieve the full spatial state vector corresponding to a reduced operation state.
		///
		/// @param reducedState State vector in the reduced operation space (e.g. |0⟩ or |1⟩).
		/// @return
		///    Full spatial state vector in the original Hilbert space, obtained by embedding the reduced state back into the full basis.
        /// @detail
		///    Reserved for visualization and spatially resolved analysis. Not used for time evolution or control dynamics,
        ///    which operate entirely within the reduced space.
        StateVector<SingleAtomFullHilbertSpace> projectToFullHilbertSpace
            (const StateVector<SingleAtomOperationHilbertSpace>& reducedState) const noexcept
        {
            return m_operationSpace.embed(*m_basisStates2D, reducedState);
		}

        /// @brief Provide the dipole matrix   
        ///
        /// @return
        ///   Dipole matrix (radial part only, assuming no variations in laser polarizationm
        ///   hence no access to states with different 'm' quantum numbers)
        static const matrix_t<ConfigType::LevelCount>& getDipoleMatrix() noexcept
        {
            return m_dipoleMatrix;
        }

        /// @brief Calculate and provide the eigen energies of all states of the operational space
        ///
        /// @return
        ///   array of Hartree energies
        static const std::array<real_t, ConfigType::LevelCount>& getHartreeEnergies() noexcept
        {
            return m_hartreeEnergies;
		}

        /// @brief Fully construct the neutral atom manifold.
        ///
        /// @details
        /// Construction proceeds in three phases:
        ///   1. Build the dipole transition matrix
        ///   2. Construct and orthonormalize full spatial basis states
        ///   3. Build the reduced energy space used for time evolution
		///   4. Calculate Hartree energies for all states
        NeutralAtomManifold()
        {
            buildDipleMatrix();
            buildFullBasisSet();
            //buildOperationSpace();
			calculateHartreeEnergies();
        }
    };
}
