// Fft.cpp
#include "Fft.h"

namespace Fft {

int PowOfTwoDataSize(int size) {
  return std::exp2(std::ceil(std::log2(size)));
}

} // namespace Fft
