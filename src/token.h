#ifndef TOKEN_H
#define TOKEN_H

namespace dstm{

struct Token {
  int time_id_;
  int doc_id_;
  int word_id_;
public:
  Token(int t,int d, int w) : time_id_(t),doc_id_(d), word_id_(w) {}
  virtual ~Token() {}
};

} // namespace lda

#endif  // TOKEN_H
