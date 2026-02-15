#pragma once
#include "core_types.h"
#include "hilbert_space/state_vector.h"

namespace KetCat
{
     /// @brief Functor to generate a 1D Quantum Harmonic Oscillator coherent state.
     ///
     /// A coherent state |α⟩ is a minimum-uncertainty wave packet that behaves
     /// most classically among quantum states. In position representation, the
     /// coherent state wavefunction is a Gaussian envelope modulated by a
     /// plane-wave phase factor:
     ///
     /// ψ(x) = exp( - (x - x₀)² / (2 σ²) ) · exp( i p₀ x / ℏ )
     ///
     /// where:
     /// - x₀ is the initial spatial displacement (center of the packet)
     /// - p₀ is the initial momentum
     /// - σ = √(ℏ / (2 m ω)) is the ground-state width of the harmonic oscillator
     /// - ℏ is the reduced Planck constant
     ///
     /// The resulting state has:
     /// - ⟨x⟩ = x₀
     /// - ⟨p⟩ = p₀
     /// - minimal Heisenberg uncertainty (Δx Δp = ℏ / 2)
     ///
     /// This implementation discretizes the wavefunction on a 1D spatial grid
     /// and returns a normalized StateVector.
     ///
     /// @tparam Dim  Number of discrete spatial grid points
    template<dimension_t Dim>
	struct CoherentStateGaussian
    {
        /// @brief Generate a discretized coherent state wavefunction.
        ///
        /// Constructs a Gaussian wave packet centered at position x₀ with
        /// momentum p₀, sampled on a uniform spatial grid with spacing dx.
        ///
        /// The wavefunction is normalized numerically using the grid spacing.
        ///
		/// @param hBar   Reduced Planck constant ℏ
        /// @param x0     Initial position ⟨x⟩ (center of the wave packet)
        /// @param p0     Initial momentum ⟨p⟩
        /// @param dx     Spatial grid spacing
        /// @param m      Particle mass (default = 1)
        /// @param omega  Harmonic oscillator angular frequency (default = 1)
		/// @param sigmaOverride Optional override for the Gaussian width σ (default = -1.0, which uses the standard σ, which is √(ℏ / (2 m ω)))
        ///
        /// @return Normalized discrete coherent state |ψ⟩
        constexpr StateVector<InfiniteHilbertSpace<Dim>> operator()(
            real_t hBar,
            real_t x0,        
            real_t p0,       
            real_t dx,
            real_t m, 
            real_t omega,
			real_t sigmaOverride = -1.0
        ) const noexcept
        {
            /// Width of the ground-state Gaussian:
            /// σ = √(ℏ / (2 m ω))
			/// In case of an override, use that value instead.
            const real_t sigma =
                (sigmaOverride == -1.0 ?
                    ConstexprMath::sqrt(hBar / (2.0 * m * omega))
                    : sigmaOverride
				);

            StateVector<InfiniteHilbertSpace<Dim>> Psi{};

            for (dimension_t i = 0; i < Dim; ++i)
            {
                const real_t x = (i + 1) * dx;

                /// Gaussian envelope centered at x₀
                /// exp( - (x - x₀)² / (2 σ²) )
                const real_t Gauss =
                    ConstexprMath::exp<20>(
                        -(x - x0) * (x - x0)
                        / (2.0 * sigma * sigma)
                    );

                /// Plane-wave phase factor argument: p₀ x / ℏ
                const real_t Phase = p0 * x / hBar;

                /// Construct complex wavefunction value:
                /// ψ(x) = gauss · (cos(phase) + i sin(phase))
                Psi[i] = cplx_t(
                    Gauss * ConstexprMath::cos(Phase),
                    Gauss * ConstexprMath::sin(Phase)
                );
            }

            Psi.normalize(dx);

            return Psi;
        }
    };
}
