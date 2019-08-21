#ifndef DSTM_H
#define DSTM_H

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <ctime>

#include "token.h"
#include "random.h"

namespace dstm {

class DSTM {
  public:
    DSTM(int* doc_num, int super_topic_num,int sub_topic_num, int word_num,int time_num,
      std::vector<Token>& tok_train_lis,std::vector<Token>& tok_test_lis) {
    D_ = doc_num;
    S_ = super_topic_num;
    K_ = sub_topic_num;
    W_ = word_num;
    T_ = time_num;
    tokens_ = tok_train_lis;
    tokens_test = tok_test_lis;

    ptr_random_ = new Random();

    Init();
  }

  double*** phi;
  double*** beta_;
  double*** alpha2_;
  int T_;                         // # of time
  int* D_;                         // # of document
  int S_;
  int K_;                       // # of topic


  // hyper parameter
  double** alpha1_;

  double**** theta_sk;
  double*** theta_s;


  std::vector<Token> tokens_;
  std::vector<Token> tokens_test;
  int* z_super;
  int* z_sub;




  ~DSTM() {
    try {
      for (int t = 0; t < T_; t++) {
        for (int k=0 ; k<K_; k++){
          delete [] ntkv[t][k];
          delete [] beta_[t][k];
        }
        delete [] ntkv[t];
        delete [] ntk[t];
        delete [] beta_[t];

        for (int d=0;d<D_[t];d++){
          for (int s=0;s<S_;s++){
              delete [] ntdsk[t][d][s];
          }
          delete [] ntds[t][d];
          delete [] ntdsk[t][d];
        }
        delete [] ntds[t];
        delete [] ntdsk[t];
        delete [] ntd[t];
      }

      delete [] ntkv;
      delete [] ntk;
      delete [] ntds;
      delete [] ntdsk;
      delete [] ntd;
      delete [] beta_;

      delete [] z_super;
      delete [] z_sub;
      delete [] p_;

      delete ptr_random_;


    for (int s = 0; s<S_; s++){
        for (int k = 0; k<K_; k++){
            delete [] histogram_ndsk[s][k];
            }
        delete [] histogram_nds[s];
        delete [] histogram_ndsk[s];
        delete [] nonZeroLimits[s];

        }
    delete [] histogram_nd;    
    delete [] histogram_nds;
    delete [] histogram_ndsk;
    delete [] nonZeroLimits;
      }

    catch (...) {
      std::cerr << "~DSTM(): Out of memory" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  void Update(int t);
  void Update_phi(int t);
  void Update_Theta_s(int t);
  void Update_Theta_sk(int t);
  void Update_alpha2(int t);
  void Update_alpha1(int t);
  void Update_beta(int t);

  double PPL(int t);

 private:
  int W_;                        // # of unique word

  int*** ntkv ;
  int*** ntds;
  int**** ntdsk;

  int** ntk;
  int** ntd;

  double* Btk;

  int*** histogram_ndsk;
  int** histogram_nds;
  int* histogram_nd;

  int* begin_train;
  int* end_train;

  int* begin_test;
  int* end_test;

  int** nonZeroLimits;

  double* p_;

  Random* ptr_random_;

  void Init();
  void _Init();

  int SelectNextTopic(const Token t);

  void Resample(std::size_t token_id);

};

} // namespace
#endif
