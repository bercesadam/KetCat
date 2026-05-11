#pragma once

namespace KetCat::Visu
{
	/// @brief Enum to specify whether to use phase encoding in visualization
	enum class UsePhaseEncoding
	{
		YES,
		NO
	};

	/// @brief Enum to specify whether to clear the screen before updating visualization
	enum class ClearScreen
	{
		YES,
		NO
	};

	/// @brief Enum to enable visualization of real and imaginary parts
	enum class ShowComplexParts
	{
		YES,
		NO
	};

	/// @brief Enum to specify whether to show potential in visualization
	enum class ShowPotential
	{
		YES,
		NO
	};

	/// @brief Check if an enum flag is enabled
	template<typename EnumType>
	constexpr inline bool enabled(EnumType e, EnumType yes = EnumType::YES) noexcept
	{
		return e == yes;
	}
}
