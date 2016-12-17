//============================================================================
// Name        : main.cpp
// Date        : 17.12.2016
// Author      : Andrey Popov
// Copyright   : All rights reserved
//============================================================================

#include "Fft.h"

int main(int argc, char* argv[]) {
  std::vector<std::complex<double>> data = {
    1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  };
  Fft::RealFastFourierTransform(data);
  return 0;
}
