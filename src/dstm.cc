#include "dstm.h"
#include <assert.h>
#include <math.h>
#include <boost/math/special_functions/digamma.hpp>

namespace dstm {

void DSTM::Init() {
    try {

    //n_tkv
    ntkv = new int**[T_];
    for (int t=0; t<T_; ++t){
        ntkv[t] = new int*[K_];
        for (int k = 0; k < K_; ++k) {
            ntkv[t][k] = new int[W_];
            for (int v=0; v<W_; ++v){
                ntkv[t][k][v] = 0;
            }
        }
    }

    //phi
    phi = new double**[T_];
    for (int t=0;t<T_; ++t){
        phi[t] = new double*[K_];
        for (int k = 0; k < K_; ++k){
            phi[t][k] = new double[W_];
            for (int v=0; v<W_; ++v){
                phi[t][k][v] = 0;
            }
        }
    }

    //theta
    theta_sk = new double***[T_];
    theta_s = new double**[T_];
    for (int t=0;t<T_; t++){
        theta_sk[t] = new double**[D_[t]];
        theta_s[t] = new double*[D_[t]];
        for (int d = 0; d < D_[t]; ++d) {
            theta_sk[t][d] = new double*[S_];
            theta_s[t][d] = new double[S_];
            for (int s = 0; s <  S_; ++s){
                theta_sk[t][d][s] = new double[K_];
                theta_s[t][d][s] = 0.0;
                for (int k = 0; k < K_; ++k) {
                      theta_sk[t][d][s][k] = 0.0;
                }
            }
        }
    }

    //n_tk
    ntk = new int*[T_];
    for (int t = 0; t < T_; t++){
        ntk[t]=new int[K_];
        for (int k = 0; k < K_; ++k){
            ntk[t][k] = 0;
        }
    }

    //n_tds
    ntds = new int**[T_];
    for (int t=0;t<T_; t++){
        ntds[t] = new int*[D_[t]];
        for (int d = 0; d < D_[t]; ++d) {
            ntds[t][d] = new int[S_];
            for (int s=0; s<S_;++s){
                ntds[t][d][s] = 0;
            }
        }
    }

    //ntdsk
    ntdsk = new int***[T_];
    for (int t=0;t<T_; t++){
        ntdsk[t] = new int**[D_[t]];
        for (int d = 0; d < D_[t]; ++d) {
            ntdsk[t][d] = new int*[S_];
            for (int s=0; s<S_;++s){
                ntdsk[t][d][s] = new int[K_];
                for (int k=0;k<K_;++k){
                    ntdsk[t][d][s][k] = 0;
                }
            }
        }
    }

    //ntd
    ntd = new int*[T_];
    for (int t = 0; t < T_; t++){
        ntd[t]=new int[D_[t]];
        for (int d = 0; d < D_[t]; ++d){
            ntd[t][d] = 0;
        }
    }

    //z ::
    z_super = new int[tokens_.size()];
    z_sub  = new int[tokens_.size()];
    for (std::size_t i = 0; i<tokens_.size(); ++i){
        z_super[i] = 0;
        z_sub[i] = 0;
    }

    //p ::
    p_ = new double[S_*K_];
    for (int i = 0; i < S_*K_; ++i)
        p_[i] = 0.0;
    
    //beta ::
    beta_ = new double**[T_];
    for (int t=0;t<T_; t++){
        beta_[t] = new double*[K_];
        for (int k = 0; k < K_; ++k) {
            beta_[t][k] = new double[K_];
            for (int k_=0; k_ < K_; ++k_){
                if (k==k_){
                    beta_[t][k][k_] = 100.0;
                }
                else{
                    beta_[t][k][k_] = 0.1;
                }
            }
        }
    }
    
    //alpha
    alpha1_ = new double* [T_];
    alpha2_ = new double**[T_];
    for(int t=0; t<T_; t++){
        alpha1_[t] = new double[S_];
        alpha2_[t] = new double*[S_];
        for (int s=0; s<S_; ++s){
            alpha1_[t][s] = 0.1;
            alpha2_[t][s] = new double[K_];
            for (int k=0;k<K_;++k){
                alpha2_[t][s][k] = 1.0;
            }
        }
    }

    int maxTokens = 11464;

    histogram_nd = new int[maxTokens];

    for (int i = 0; i<maxTokens; ++i){
        histogram_nd[i] = 0;
    }

    histogram_nds = new int*[S_];
    histogram_ndsk = new int**[S_];

    for (int s = 0; s<S_; ++s){
        //make nds histgram
        histogram_nds[s] = new int[maxTokens];
        histogram_ndsk[s] = new int*[K_];
        for (int k = 0; k<K_; ++k){
            //ndsk ndsk histgram
            histogram_ndsk[s][k] = new int[maxTokens];
            for (int i = 0; i<maxTokens; ++i){
                histogram_nds[s][i] = 0;
                histogram_ndsk[s][k][i] = 0;
            }
        }
    }

    nonZeroLimits = new int*[S_];
    for (int s=0; s<S_;++s){
        nonZeroLimits[s] = new int[K_];
        for (int k=0; k<K_;++k){
            nonZeroLimits[s][k] = -1;
        }
    }

    Btk = new double[W_];


    begin_train = new int[T_];
    end_train = new int[T_];

    begin_test = new int[T_];
    end_test = new int[T_];

    for (int t=0;t<T_; ++t){
        begin_train[t] = -1;
        begin_test[t] = -1;
    }

    for (std::size_t i = 0; i < tokens_.size(); ++i){
        Token tok = tokens_[i];
        int t = tok.time_id_;
        if (begin_train[t]==-1){
            begin_train[t] = i;
        }
        end_train[t] = i;
    }

    for (std::size_t i = 0; i < tokens_test.size(); ++i){
        Token tok = tokens_test[i];
        int t = tok.time_id_;
        if (begin_test[t]==-1){
            begin_test[t] = i;
        }
        end_test[t] = i;
    }
    }


    catch (...) {
        std::cerr << "Init(): Out of memory" << std::endl;
        exit(EXIT_FAILURE);
    }

    _Init();
}

void DSTM::_Init() {
    for (std::size_t i = 0; i<tokens_.size(); ++i) {
        Token tok = tokens_[i];

        int t = tok.time_id_;
        int d = tok.doc_id_;
        int v = tok.word_id_;

        std::size_t assign_s = static_cast<std::size_t>((*ptr_random_).gen(S_));
        std::size_t assign_k = static_cast<std::size_t>((*ptr_random_).gen(K_));

        ntkv[t][assign_k][v]++;
        ntds[t][d][assign_s]++;
        ntdsk[t][d][assign_s][assign_k]++;

        ntk[t][assign_k]++;
        ntd[t][d]++;

        z_super[i] = assign_s;
        z_sub[i] = assign_k;
    }
}

int DSTM::SelectNextTopic(Token tok) {
    int t = tok.time_id_;
    int d = tok.doc_id_;
    int v = tok.word_id_;

    double alpha2_sum[S_];

    for(int s=0; s <S_; ++s){
        alpha2_sum[s] = 0.0;
        for (int k=0; k < K_; ++k){
            alpha2_sum[s] += alpha2_[t][s][k];
        }
    }


    if (t>0){
        for (int s=0;s<S_;++s){
            for (int k = 0; k < K_; ++k) {
                int tmp = s*K_+k;
                double tmp1 = 0.0;
                double tmp2 = 0.0;

                for (int k_=0; k_<K_; ++k_){
                    tmp1 += beta_[t][k][k_]*phi[t-1][k_][v];
                    tmp2 += beta_[t][k][k_];
                }

                p_[tmp] = (ntkv[t][k][v] + tmp1)/(ntk[t][k] + tmp2)
                       *(ntds[t][d][s] + alpha1_[t][s])
                       *(ntdsk[t][d][s][k] + alpha2_[t][s][k])/(ntds[t][d][s] + alpha2_sum[s]);

                if (tmp != 0) p_[tmp] += p_[tmp - 1];
            }
        }
    }
    else{
        for (int s=0; s<S_; ++s){
            for (int k = 0; k < K_; ++k) {
                int tmp = s*K_+k;
                p_[tmp] = (ntkv[t][k][v] + 0.1)/(ntk[t][k] + W_*0.1)
                            *(ntds[t][d][s] + alpha1_[t][s])
                            *(ntdsk[t][d][s][k] + alpha2_[t][s][k])/(ntds[t][d][s] + alpha2_sum[s]);

                if (tmp != 0) p_[tmp] += p_[tmp - 1];
            }
        }
    }

    int max = S_*K_;
    double u = (*ptr_random_).gen(1.0) * p_[max-1];

    for (int i = 0; i < max; ++i) {
        if (u < p_[i]) return i;
    }
    return S_*K_ - 1;
}

inline void DSTM::Resample(std::size_t token_id) {
    Token tok = tokens_[token_id];
    int assign_s = z_super[token_id];
    int assign_k = z_sub[token_id];

    int t = tok.time_id_;
    int d = tok.doc_id_;
    int v = tok.word_id_;

    ntkv[t][assign_k][v]--;
    ntds[t][d][assign_s]--;
    ntdsk[t][d][assign_s][assign_k]--;
    ntk[t][assign_k]--;

    int assign = SelectNextTopic(tok);
    assign_s = assign/K_;
    assign_k = assign%K_;

    ntkv[t][assign_k][v]++;
    ntds[t][d][assign_s]++;
    ntdsk[t][d][assign_s][assign_k]++;
    ntk[t][assign_k]++;

    z_super[token_id] = assign_s;
    z_sub[token_id] = assign_k;
}

void DSTM::Update(int t) {
    for (std::size_t i = begin_train[t]; i<=end_train[t]; ++i){
        DSTM::Resample(i);
    }
    DSTM::Update_phi(t);
    DSTM::Update_Theta_s(t);
    DSTM::Update_Theta_sk(t);
}

void DSTM::Update_phi(int t){
    if (t>0){
        for (int k = 0; k < K_; ++k) {
            double sum = 0.0;
            for (int v = 0; v < W_; ++v) {
                double tmp = 0.0;
                for (int k_ = 0; k_ < K_; ++k_) {
                    tmp += beta_[t][k][k_]*phi[t-1][k_][v];
                }
                phi[t][k][v] = tmp + ntkv[t][k][v];
                sum += phi[t][k][v];
            }
            // normalize
            double sinv = 1.0 / sum;
            for (int v = 0; v < W_; ++v)
                phi[t][k][v] *= sinv;
        }
    }
    else{
        for (int k = 0; k < K_; ++k) {
            double sum = 0.0;
            for (int v = 0; v < W_; ++v) {
                phi[t][k][v] = 0.1 + ntkv[t][k][v];
                sum += phi[t][k][v];
            }
            // normalize
             double sinv = 1.0 / sum;
            for (int v = 0; v < W_; ++v)
                phi[t][k][v] *= sinv;
        }
    }
}
void DSTM::Update_alpha1(int t) {
    int maxTokens = 11464; //pending

    for (int i = 0; i<maxTokens; ++i){
        histogram_nd[i] = 0;
    }

    for (int s = 0; s<S_; ++s){
        for (int i = 0; i<maxTokens; ++i){
            histogram_nds[s][i] = 0;
        }
    }

    for (int d=0; d<D_[t];++d){
        histogram_nd[ntd[t][d]] += 1;
    }

    for (int d=0; d<D_[t];++d){
        for (int s=0; s<S_;++s){
            histogram_nds[s][ntds[t][d][s]] += 1;
        }
    }

    //iteration part
    for (int i=0 ; i<20; ++i){
        double alpha1_sum = 0.0;
        for(int s=0; s <S_; ++s){
            alpha1_sum += alpha1_[t][s];
        }

        double deli = 0.0;
        double currentDigamma_deli = 0.0;
        for (int i=1; i<maxTokens; ++i){
            currentDigamma_deli += 1.0/(alpha1_sum+i-1);
            deli += histogram_nd[i]*currentDigamma_deli;
        }

        for (int s=0; s<S_; ++s){
            double nume = 0.0;
            double currentDigamma_nume = 0.0;

            for (int i=1; i<maxTokens; ++i){
                currentDigamma_nume += 1.0/(alpha1_[t][s]+i-1);
                nume += histogram_nds[s][i]*currentDigamma_nume;
            }
          alpha1_[t][s] *= nume/deli;
        }
    }
}

void DSTM::Update_alpha2(int t) {
    int maxTokens = 11464; //pending

    for (int s = 0; s<S_; ++s){
        //make nds histgram
        for (int k = 0; k<K_; ++k){
            //make ndsk histgram
            for (int i = 0; i<maxTokens; ++i){
                histogram_nds[s][i] = 0;
                histogram_ndsk[s][k][i] = 0;
            }
        }
    }

    for (int d=0; d<D_[t];++d){
        for (int s=0; s<S_;++s){
            histogram_nds[s][ntds[t][d][s]] += 1;
            for (int k=0; k<K_;++k){
                histogram_ndsk[s][k][ntdsk[t][d][s][k]] += 1;
            }
        }
    }

    for (int s=0; s<S_;++s){
        for (int k=0; k<K_;++k){
            nonZeroLimits[s][k] = -1;
        }
    }

    for (int s=0; s<S_; ++s) {
        for (int k=0; k<K_; ++k) {
            for (int i = 1; i < maxTokens; ++i) {
                if (histogram_ndsk[s][k][i] > 0){
                    nonZeroLimits[s][k] = i;
                }
            }
        }
    }

    //iteration part
    for (int i=0 ; i<20; ++i){
        double alpha2_sum[S_];
        for(int s=0; s <S_; ++s){
            alpha2_sum[s] = 0.0;
            for (int k=0; k < K_; ++k){
                alpha2_sum[s] += alpha2_[t][s][k];
            }
        }

        for (int s=0; s<S_; ++s){
            double deli = 0.0;
            double currentDigamma_deli = 0.0;

            for (int i=1; i<maxTokens; ++i){
                currentDigamma_deli += 1.0/(alpha2_sum[s]+i-1);
                deli += histogram_nds[s][i]*currentDigamma_deli;
            }
            
            for (int k=0; k<K_;++k){

                if (nonZeroLimits[s][k]==-1){
                    alpha2_[t][s][k] = 0.000001;
                    continue;
                }

                double nume = 0.0;
                double currentDigamma_nume = 0.0;

                for (int i=1; i <= nonZeroLimits[s][k]; ++i){
                    currentDigamma_nume += 1.0/(alpha2_[t][s][k]+i-1);
                    nume += histogram_ndsk[s][k][i]*currentDigamma_nume;
                }
               
                alpha2_[t][s][k] *= nume/deli;
            }
        }
    }
}

void DSTM::Update_beta(int t){
    for (int i = 0 ; i<20; ++i){
        for (int k=0; k<K_;++k){
            double betak_sum = 0.0;

            for (int k_=0; k_<K_;++k_){
                betak_sum += beta_[t][k][k_];
            }

            double deli = boost::math::digamma(ntk[t][k]+betak_sum)
                         -boost::math::digamma(betak_sum);

            for (int v=0; v<W_;++v){
                double beta_tmp = 0.0;
                for (int k_=0; k_<K_;++k_){
                    beta_tmp += beta_[t][k][k_]*phi[t-1][k_][v];
                }
                Btk[v] =  boost::math::digamma(ntkv[t][k][v]+beta_tmp)
                          -boost::math::digamma(beta_tmp);
            }

            for (int k_=0; k_<K_;++k_){
                double nume =0.0;
                for (int v=0; v<W_;++v){
                    nume += phi[t-1][k_][v]*Btk[v];
                }
                beta_[t][k][k_] = nume/deli*beta_[t][k][k_];
            }
        }
    }
}

void DSTM::Update_Theta_s(int t) {
    for (int d = 0; d < D_[t]; ++d) {
        double sum = 0.0;
        for (int s = 0; s < S_; ++s) {
            theta_s[t][d][s] = alpha1_[t][s] + ntds[t][d][s];
            sum += theta_s[t][d][s];
        }
        // normalize
        double sinv = 1.0 / sum;
        for (int s = 0; s < S_; ++s) {
            theta_s[t][d][s] *= sinv;
        }
    }
}

void DSTM::Update_Theta_sk(int t) {
    for (int s = 0; s < S_; ++s) {
        for (int d = 0; d < D_[t]; ++d) {
            double sum = 0.0;
            for (int k = 0; k < K_; ++k) {
                theta_sk[t][d][s][k] = alpha2_[t][s][k] + ntdsk[t][d][s][k];
                sum += theta_sk[t][d][s][k];
            }
            // normalize
            double sinv = 1.0 / sum;
            for (int k = 0; k < K_; ++k) {
                theta_sk[t][d][s][k] *= sinv;
            }
        }
    }
}

double DSTM::PPL(int t) {
    double log_per = 0.0;
    double nt = double(end_test[t]-begin_test[t]);

    for (std::size_t i = begin_test[t]; i<=end_test[t]; ++i){
        Token tok = tokens_test[i];
        int d = tok.doc_id_;
        int v = tok.word_id_;

        double tmp = 0.0;
        for(int s = 0; s < S_; ++s){
            for (int k=0; k<K_; ++k){
                tmp += theta_s[t][d][s]*theta_sk[t][d][s][k]*phi[t][k][v];
            }
        }
        log_per -= log(tmp);
    }
    return exp(log_per/nt);
    }
} 
// namespace dstm
