//============================================================================
// Name        : Fft.h
// Date        : 17.12.2016
// Author      : Andrey Popov
// Copyright   : All rights reserved
//============================================================================

#include "Fft.h"

namespace Fft {

int PowOfTwoDataSize(int size) {
  return std::exp2(std::ceil(std::log2(size)));
}

} // namespace Fft
