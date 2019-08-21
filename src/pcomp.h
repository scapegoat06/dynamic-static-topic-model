#ifndef PCOMP_H
#define PCOMP_H

namespace dstm {

struct Pcomp {
  int id;
  double prob;
};

class LessProb {
public:
  bool operator()(const Pcomp& a, const Pcomp& b) const {
    return b.prob < a.prob;
  }
};

} // namespace lda

#endif  // PCOMP_H
