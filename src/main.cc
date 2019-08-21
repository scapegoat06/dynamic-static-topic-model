#include <cstring>
#include "dstm.h"
#include "pcomp.h"
#include <math.h>
#include <sys/stat.h>
#include <stdio.h>
#include <unistd.h>


namespace dstm {

    inline int ParseLine(const std::string& line, std::vector<Token>& tokens) {
        std::istringstream is(line);
        int time_id =0;
        int doc_id = 0;
        int word_id = 0;
        int count = 0;
        is >> time_id >>doc_id >> word_id >> count;

        if (!doc_id || !word_id || !count) {
            std::cerr << "parse error";
            return -1;
        }

        for (int i = 0; i < count; ++i) {
            tokens.push_back(Token(time_id -1 ,doc_id - 1, word_id - 1));
        }

        return 0;
    }

    int ReadBOWData(const std::string& file, std::vector<Token>& tokens,
                    int& T,int D[], int& W, int& N) {

        std::istream *ifs;

        if (file == "-") {
            ifs = &std::cin;
        }
        else {
            ifs = new std::ifstream(file.c_str());
        }

        if (!*ifs) {
            std::cerr << "Cannot open: " << file << std::endl;
            return -1;
        }

        std::size_t line_num = 0;
        std::string line;

        // Get feature size
        std::getline(*ifs, line);
        T = atoi(line.c_str());
        std::getline(*ifs,line);
        std::string token;
        std::istringstream stream(line);

        int t = 0;
        while(std::getline(stream,token,',')){
            int temp=atoi(token.c_str()); //stof(string str) : string to float
            D[t] = temp;
            t += 1;
        }

        std::getline(*ifs, line);
        W = atoi(line.c_str());
        std::getline(*ifs, line);
        N = atoi(line.c_str());

        line_num += 4;

        if (N <= 0) {
            std::cerr << "Invalid # of N" << std::endl;
            return -1;
        }

        for (int i = 0; i < N; ++i) {
            std::getline(*ifs, line);
            if (line[0] == '#') continue;
            if (ParseLine(line, tokens) == - 1) {
                std::cerr << " line: " << line_num;
                return -1;
            }
            line_num++;
        }
        if (file != "-") delete ifs;

        return 0;
    }

    int ReadVocabData(const std::string& file,
                    std::vector<const char*>& words, int& W) {
        std::istream *ifs;

        if (file == "-") {
            ifs = &std::cin;
        }
        else {
            ifs = new std::ifstream(file.c_str());
        }

        if (!*ifs) {
            std::cerr << "Cannot open: " << file << std::endl;
            return -1;
        }

        std::string line;
        for (int i = 0; i < W; ++i) {
            std::getline(*ifs, line);
            if (line[0] == '#') continue;
            char *tmp = new char[line.size()+1];
            std::strcpy(tmp, line.c_str());
            words[i] = tmp;
        }
        if (file != "-") delete ifs;

        return 0;
    }

    void PrintWordTopic(double** phi, int K, int W,
                        std::vector<const char *>& words) {
        for (int k = 0; k < K; ++k) {
            std::cout << "topic: " << k << std::endl;
            std::vector<Pcomp> ps(W);

            for (int w = 0; w < W; ++w) {
                Pcomp pc;
                pc.id = w;
                pc.prob = phi[k][w];
                ps[w] = pc;
            }

            std::sort(ps.begin(), ps.end(), LessProb());

            // print top 10 words
            for (int i = 0; i < 10; ++i) {
                Pcomp p = ps[i];
                std::cout << words[p.id] << ' ' << p.prob << std::endl;
            }
            std::cout << std::endl;
        }
    }

    void save_phi(double*** phi,int T,int K,int V,std::string fname){
        std::string filename = fname;
        std::ofstream writing_file;
        writing_file.open(filename, std::ios::out);

        writing_file << T << std::endl;
        writing_file << K << std::endl;
        writing_file << V << std::endl;

        for (int t = 0; t<T; ++t){
            for (int k = 0; k < K; ++k) {
                for (int v = 0; v < V; ++v) {
                    writing_file << t << "," << k << "," << v << "," << phi[t][k][v]<< std::endl;
                }
            }
        }
    }

    void save_theta_s(double*** theta_s,int T,int* D,int S,std::string fname){
        std::string filename = fname;
        std::ofstream writing_file;
        writing_file.open(filename, std::ios::out);

        writing_file << T << std::endl;
        for (int t = 0; t<T; ++t){
            if (t==0){
                writing_file << D[t];
            }
            else{
                writing_file << "," <<D[t];
            }
        }
        writing_file << std::endl;
        writing_file << S << std::endl;

        for (int t = 0; t<T; ++t){
            for (int d = 0; d < D[t]; ++d) {
                for (int s = 0; s < S; ++s) {
                    writing_file << t << "," << d << "," << s << "," << theta_s[t][d][s]<< std::endl;
                }
            }
        }
    }

    void save_theta_sk(double**** theta_sk,int T,int* D,int S,int K,std::string fname){
            std::string filename = fname;
            std::ofstream writing_file;
            writing_file.open(filename, std::ios::out);

            writing_file << T << std::endl;
            for (int t = 0; t<T; ++t){
                if (t==0){
                    writing_file << D[t];
                }
                else{
                    writing_file << "," <<D[t];
                }
            }
            writing_file << std::endl;
            writing_file << S << std::endl;
            writing_file << K << std::endl;

            for (int t = 0; t<T; ++t){
                for (int d = 0; d < D[t]; ++d) {
                    for (int s = 0; s < S; ++s) {
                        for (int k=0; k < K; ++k){
                            writing_file << t << "," << d << "," << s << "," <<k<<","<< theta_sk[t][d][s][k]<< std::endl;
                        }
                    }
                }
            }
        }

    void save_beta(double*** beta,int T,int K,std::string fname){
        std::string filename = fname;
        std::ofstream writing_file;
        writing_file.open(filename, std::ios::out);

        writing_file << T << std::endl;
        writing_file << K << std::endl;

        for (int t = 0; t<T; ++t){
            for (int k = 0; k < K; ++k) {
                for (int k_ = 0; k_ < K; ++k_) {
                    writing_file << t << "," << k << "," << k_ << "," << beta[t][k][k_]<< std::endl;
                }
            }
        }
    }

    void save_alpha1(double** alpha1,int T,int S,std::string fname){
        std::string filename = fname;
        std::ofstream writing_file;
        writing_file.open(filename, std::ios::out);

        writing_file << T << std::endl;
        writing_file << S << std::endl;

        for (int t = 0; t<T; ++t){
            for (int s = 0; s < S; ++s) {
                writing_file << t << "," << s << ","  << alpha1[t][s]<< std::endl;
            }
        }
    }

    void save_alpha2(double*** alpha2,int T,int S,int K,std::string fname){
        std::string filename = fname;
        std::ofstream writing_file;
        writing_file.open(filename, std::ios::out);

        writing_file << T << std::endl;
        writing_file << S << std::endl;
        writing_file << K << std::endl;

        for (int t = 0; t<T; ++t){
            for (int s = 0; s < S; ++s) {
                for (int k = 0; k < K; ++k) {
                    writing_file << t << "," << s << "," << k << "," << alpha2[t][s][k]<< std::endl;
                }
            }
        }
    }

    void save_PPL(double* PPL,int t,int num_iter,std::string fname){
        std::string filename = fname;
        std::ofstream writing_file;
        writing_file.open(filename, std::ios::app);

        for (int i = 0; i<num_iter; ++i){
            writing_file << t << "," << i << "," << PPL[i]<< std::endl;
        }
    }

}//end namespace dstm


std::string IntToString(int number){
  std::stringstream ss;
    ss << number;
      return ss.str();
}


int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " bow_file vocab_file" << std::endl;
        return -1;
    }

    std::string bow_train = argv[1];
    std::string bow_test = argv[2];
    std::string vocab_file = argv[3];

    const int S = atoi(argv[4]);
    const int K = atoi(argv[5]);
    const int num_iter = atoi(argv[6]);

    std::vector<dstm::Token> tokens_train;
    std::vector<dstm::Token> tokens_test;

    int D[210]; //pending
    int W = 0;
    int N = 0;
    int T = 0;

    // Read Bow-train file
    if (dstm::ReadBOWData(bow_train, tokens_train,T, D, W, N) != 0) {
        std::cerr << "Cannot read" << std::endl;
        return -1;
    }

    // Read Bow-test file
    if (dstm::ReadBOWData(bow_test, tokens_test,T, D, W, N) != 0) {
        std::cerr << "Cannot read" << std::endl;
        return -1;
    }

    // Read vocabrary file
    std::vector<const char*> words(W);
    if (dstm::ReadVocabData(vocab_file, words, W) != 0) {
        std::cerr << "Cannot read" << std::endl;
        return -1;
    }

    //Data details
    std::cout <<"Number of TimeStep : "<< T <<std::endl;
    for (int t=0;t<T;++t){
        std::cout <<"Number of Documents : "<< D[t] <<std::endl;
    }
    std::cout <<"Vocablary  Size : "<< W <<std::endl;
    std::cout <<"Total Token Size : "<< tokens_train.size() <<std::endl;


    std::cout << "Start " << std::endl;

    //Get current time
    time_t now = time(NULL);
    struct tm *pnow = localtime(&now);

    std::string year = IntToString(pnow->tm_year + 1900);
    std::string mon = IntToString(pnow->tm_mon + 1);
    std::string day = IntToString(pnow->tm_mday);
    std::string hour = IntToString(pnow->tm_hour);
    std::string min = IntToString(pnow->tm_min);

    std::string Snum = IntToString(S);
    std::string Knum = IntToString(K);

    std::string cd = "../output/";
    std::string joint = "S"+Snum+"K"+Knum+":";
    std::string newname = cd+joint+year+"-"+mon+"-"+day+"-"+hour+"-"+min;

    mkdir(newname.c_str(), 0775);

    //Start Model Process

    dstm::DSTM dstm(D,S,K,W,T,tokens_train,tokens_test);

    for (int t = 0; t<T; ++t){
        double PPL[num_iter];
        for (int i = 0; i <= num_iter; ++i) {
            std::cout <<i<<std::endl;
            dstm.Update(t);
            if (i >15 and (i+1)%5==0){
                    //dstm.Update_alpha1(t);
                    dstm.Update_alpha2(t);
            }
            
            //PPL calculate
            double ppl = dstm.PPL(t);
            PPL[i] = ppl;
            if (i%100==0){
                    std::cout <<t<<" "<<i<<" "<< ppl << std::endl;
            }
            
            //Update beta
            if (t>0 and (i+1)>=50 and (i+1)%5==0){
                dstm.Update_beta(t);
            }
        }

        dstm::PrintWordTopic(dstm.phi[t], K, W, words);

        double min_ppl = PPL[0];
        int min_index = 0;
        for (int i = 0; i < num_iter; ++i) {
                if (PPL[i]<min_ppl){
                        min_ppl = PPL[i];
                        min_index = i;
                        }
                }

        std::cout <<"t="<<t<<" "<<"i="<<min_index<<" "<<"min_ppl:"<<min_ppl<<std::endl;
        dstm::save_PPL(PPL,t,num_iter,(newname+"/"+"ppl.txt").c_str());
    }

    dstm::save_phi(dstm.phi,T,K,W,(newname+"/"+"phi.txt").c_str());
    dstm::save_beta(dstm.beta_,T,K,(newname+"/"+"beta.txt").c_str());
    dstm::save_alpha1(dstm.alpha1_,T,S,(newname+"/"+"alpha1_.txt").c_str());
    dstm::save_alpha2(dstm.alpha2_,T,S,K,(newname+"/"+"alpha2_.txt").c_str());
    dstm::save_theta_s(dstm.theta_s,T,D,S,(newname+"/"+"theta_s.txt").c_str());
    dstm::save_theta_sk(dstm.theta_sk,T,D,S,K,(newname+"/"+"theta_sk.txt").c_str());

    for (std::size_t i = 0; i < words.size(); ++i)
        delete [] words[i];

    return 0;
}
