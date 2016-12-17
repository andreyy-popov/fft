//============================================================================
// Name        : Fft.h
// Date        : 17.12.2016
// Author      : Andrey Popov
// Copyright   : All rights reserved
//============================================================================

#ifndef _FFT_H_
#define _FFT_H_

#include <list>
#include <vector>
#include <complex>
#include <algorithm>
#include <functional>
#include <utility>
#include <thread>
#include <future>
#include <cmath>

namespace Fft {

template <typename T>
class Weight {
  std::vector<std::complex<T>> w;
public:
  // Конструктор.
  Weight(int n) : w(n / 2) {
    const T pi = std::acos(-1);
    const T theta = -2 * pi / n;
    T angle = 0;
    for (auto& x : w) {
      x = {std::cos(angle), std::sin(angle)};
      angle += theta;
    }
  }
  // Вычисление функции.
  std::complex<T> operator()(int k) const {
    return w[k];
  }
};

template <typename T>
class Butterfly {
  std::vector<std::complex<T>>& data;
  std::function<std::complex<T>(int)>& w;
public:
  // Конструктор.
  Butterfly(std::vector<std::complex<T>>& data, std::function<std::complex<T>(int)>& w) : data(data), w(w)
    {}
  // Вычисление функции.
  void operator()(int i, int j, int k) {
    const std::complex<T> x1 = data[i];
    const std::complex<T> x2 = data[j] * w(k);
    data[i] = x1 + x2;
    data[j] = x1 - x2;
  }
};

int PowOfTwoDataSize(int size);

template <typename T>
void PaddingData(std::vector<T>& data, int n) {
  data.resize(n, 0);
}

template <typename T>
void PaddingData(std::vector<T>& data) {
  PaddingData(data, Fft::PowOfTwoDataSize(data.size()));
}

template <typename T>
void ReverseBitOrder(std::vector<T>& data, int n, int cpu) {
  // Массив бит.
  typedef std::list<int> Bits;
  Bits bits;
  for (int bit = 1; bit != n; bit <<= 1) {
    bits.push_back(bit);
  }
  // Перевернутый массив бит.
  Bits revBits;
  std::reverse_copy(bits.begin(), bits.end(), std::back_inserter(revBits));
  // Функция переупорядочивания элементов.
  typedef typename Bits::iterator It;
  std::function<void(int, int, int, It, It, It, It)>
  reorder = [&data, &reorder] (int cpu, int index1, int index2, It it1, It end1, It it2, It end2) {
    if (it1 == end1) {
      if (index1 < index2) {
        std::swap(data[index1], data[index2]);
      }
    } else {
      int bit1 = *it1++;
      int bit2 = *it2++;
      if (cpu == 1) {
        reorder(cpu, index1, index2, it1, end1, it2, end2);
        reorder(cpu, index1 | bit1, index2 | bit2, it1, end1, it2, end2);
      } else {
        auto f = std::async(std::launch::async, reorder, cpu / 2, index1, index2, it1, end1, it2, end2);
        reorder(cpu / 2, index1 | bit1, index2 | bit2, it1, end1, it2, end2);
        f.get();
      }
    }
  };
  reorder(cpu, 0, 0, bits.begin(), bits.end(), revBits.begin(), revBits.end());
}

template <typename T>
void FastFourierTransform(std::vector<std::complex<T>>& data, std::function<std::complex<T>(int)>& w, int n, int cpu) {
  // Двоично-инверсная перестановка.
  ReverseBitOrder(data, n, cpu);
  // Функция вычисления "бабочки".
  Butterfly<T> butterfly(data, w);
  // Функция вычисления БПФ.
  std::function<void(int, int, int, int)>
  fft = [&data, &butterfly, &fft] (int cpu, int start, int end, int step) {
    if ((end - start) != 1) {
      int mid = (start + end) / 2;
      int dubstep = step * 2;
      if (cpu == 1) {
        fft(cpu, start, mid, dubstep);
        fft(cpu, mid, end, dubstep);
      } else {
        auto f = std::async(std::launch::async, fft, cpu / 2, start, mid, dubstep);
        fft(cpu / 2, mid, end, dubstep);
        f.get();
      }
      for (int i = start, j = mid, k = 0; j < end; ++i, ++j, k += step) {
        butterfly(i, j, k);
      }
    }
  };
  // Вычисляем БПФ.
  fft(cpu, 0, n, 1);
}

template <typename T>
void FastFourierTransform(std::vector<std::complex<T>>& data, int n, int cpu) {
  std::function<std::complex<T>(int)> w = Weight<T>(n);
  FastFourierTransform(data, w, n, cpu);
}

template <typename T>
void FastFourierTransform(std::vector<std::complex<T>>& data, int cpu) {
  FastFourierTransform(data, data.size(), cpu);
}

template <typename T>
void FastFourierTransform(std::vector<std::complex<T>>& data) {
  FastFourierTransform(data, data.size(), std::thread::hardware_concurrency());
}

template <typename T>
void InverseFastFourierTransform(std::vector<std::complex<T>>& data, std::function<std::complex<T>(int)>& w, int n, int cpu) {
  T c = 1.0 / n;
  for (int i = 0; i < n; ++i) {
    data[i] = c * std::conj(data[i]);
  }
  FastFourierTransform(data, w, n, cpu);
}

template <typename T>
void InverseFastFourierTransform(std::vector<std::complex<T>>& data, int n, int cpu) {
  std::function<std::complex<T>(int)> w = Weight<T>(n);
  InverseFastFourierTransform(data, w, n, cpu);
}

template <typename T>
void InverseFastFourierTransform(std::vector<std::complex<T>>& data, int cpu) {
  InverseFastFourierTransform(data, data.size(), cpu);
}

template <typename T>
void InverseFastFourierTransform(std::vector<std::complex<T>>& data) {
  InverseFastFourierTransform(data, data.size(), std::thread::hardware_concurrency());
}

template <typename T>
void RealFastFourierTransform(std::vector<std::complex<T>>& data, std::function<std::complex<T>(int)>& w, int n, int cpu) {
  const int n2 = n / 2;
  // Объединяем четные и нечетные отсчеты.
  for (int i = 0, j = 0; i < n2; ++i, j += 2) {
    data[i] = {data[j].real(), data[j + 1].real()};
  }
  // Прямое БПФ.
  std::function<std::complex<T>(int)> w2 = [&w] (int k) { return w(k + k); };
  FastFourierTransform(data, w2, n2, cpu);
  // Разделение составляющих спектра.
  const std::complex<T> x1 = data[0];
  const std::complex<T> x2 = std::conj(x1);
  data[0] = (x1 + x2) / std::complex<T>(2);
  data[n2] = (x1 - x2) / std::complex<T>(0, 2);
  for (int i = 1, j = n2 - 1, k = n2 + 1; i < n2; ++i, --j, ++k) {
    const std::complex<T>& x1 = data[i];
    const std::complex<T>& x2 = data[j];
    data[k] = (x1 - std::conj(x2)) / std::complex<T>(0, 2);
  }
  for (int i = 1, j = n2 - 1; i < j; ++i, --j) {
    const std::complex<T> x1 = data[i];
    const std::complex<T> x2 = data[j];
    data[i] = (x1 + std::conj(x2)) / std::complex<T>(2);
    data[j] = (x2 + std::conj(x1)) / std::complex<T>(2);
  }
  // Последний шаг прямого БПФ.
  Butterfly<T> butterfly(data, w);
  for (int i = 0, j = n2; i < n2; ++i, ++j) {
    butterfly(i, j, i);
  }
}

template <typename T>
void RealFastFourierTransform(std::vector<std::complex<T>>& data, int n, int cpu) {
  std::function<std::complex<T>(int)> w = Weight<T>(n);
  RealFastFourierTransform(data, w, n, cpu);
}

template <typename T>
void RealFastFourierTransform(std::vector<std::complex<T>>& data, int cpu) {
  RealFastFourierTransform(data, data.size(), cpu);
}

template <typename T>
void RealFastFourierTransform(std::vector<std::complex<T>>& data) {
  RealFastFourierTransform(data, data.size(), std::thread::hardware_concurrency());
}

template <typename T>
void RealInverseFastFourierTransform(std::vector<std::complex<T>>& data, std::function<std::complex<T>(int)>& w, int n, int cpu) {
  const int n2 = n / 2;
  // Объединяем составляющие спектра.
  std::function<std::complex<T>(const std::complex<T>&, const std::complex<T>&, int)>
  f = [&w] (const std::complex<T>& x1, const std::complex<T>& x2, int k) {
    return
      ((x1 + std::conj(x2)) / 2.0) +
      ((x1 - std::conj(x2)) / w(k) / std::complex<T>(0, 2));
  };
  data[0] = f(data[0], data[n2], 0);
  for (int i = 1, j = n2 - 1; i < j; ++i, --j) {
    const std::complex<T> x1 = data[i];
    const std::complex<T> x2 = data[j];
    data[i] = f(x1, x2, i);
    data[j] = f(x2, x1, j);
  }
  // Обратное БПФ.
  std::function<std::complex<T>(int)> w2 = [&w] (int k) { return w(k + k); };
  InverseFastFourierTransform(data, w2, n2, cpu);
  // Разделяем четные и нечетные отсчеты.
  for (int i = n2 - 1; i >= 0; --i) {
    const std::complex<T> x = data[i];
    const int j = i + i;
    data[j] = x.real();
    data[j + 1] = x.imag();
  }
}

template <typename T>
void RealInverseFastFourierTransform(std::vector<std::complex<T>>& data, int n, int cpu) {
  std::function<std::complex<T>(int)> w = Weight<T>(n);
  RealInverseFastFourierTransform(data, w, n, cpu);
}

template <typename T>
void RealInverseFastFourierTransform(std::vector<std::complex<T>>& data, int cpu) {
  RealInverseFastFourierTransform(data, data.size(), cpu);
}

template <typename T>
void RealInverseFastFourierTransform(std::vector<std::complex<T>>& data) {
  RealInverseFastFourierTransform(data, data.size(), std::thread::hardware_concurrency());
}

} // namespace Fft

#endif // _FFT_H_
