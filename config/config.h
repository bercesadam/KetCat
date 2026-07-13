#pragma once
#include "local_space/neutral_atom_config.h"


namespace KetCat::AtomConfig
{
    using namespace SpectroscopicLetters;

    NeutralAtomTypeConfig
        <
        Element::Cs,

        456, /* Spatial discretization steps count */
        750.0, /* Spatial extent in a.u. */

        0, /* Index of the logical level 0 */
        2, /* Index of the logical level 1*/
        4, /* Index of the Rydberg level */

        QuantumNumber<6, s>,  /*0*/
        QuantumNumber<6, p>,  /*1*/
        QuantumNumber<7, s>,  /*2*/
        QuantumNumber<7, p>, /*3*/
        QuantumNumber<60, s>, /*4*/
        QuantumNumber<60, p>  /*5*/
        > Cesium_6Level;
}

namespace KetCat::SimulationConfig
{
    static constexpr natural_t SimuSaveNthFrame = 1E6;
}
