#ifndef RANDOM_H
#define RANDOM_H

namespace dstm {

class Random {
public:
  Random() {
    srand(static_cast<std::size_t>(time(NULL)));
  }

  virtual ~Random() {}

  int gen(int max) {
    double tmp = static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0);
    return static_cast<int>(tmp * max);
  }

  double gen(double max) {
    double tmp = static_cast<double>(rand()) / static_cast<double>(RAND_MAX + 1.0);
    return tmp * max;
  }
};

} // namespace lda

#endif  // RANDOM_H
