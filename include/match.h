

#ifndef INCLUDE_MATCH_H_
#define INCLUDE_MATCH_H_
#include "basic.h"
#include "kmer_TRs_interval.h"
#include "parallel.h"
#include "generateSeq.h"

struct interval_space{
    pair<int, int> inverval;
    int space;
};
struct thread_match_tr {
    vector<vector<uint32_t>> p_subset;
    vector<vector<double>> ratios_subset;
    kmer_TRs_interval& ktis;
    int kmer_len;
	thread_match_tr(kmer_TRs_interval& ktis_ref,
        vector<vector<uint32_t>> p_subset_arg,
        vector<vector<double>> ratios_subset_arg,
        int kmer_len_arg) : ktis(ktis_ref),
          p_subset(p_subset_arg),
          ratios_subset(ratios_subset_arg),
          kmer_len(kmer_len_arg) {}

    // 显式删除拷贝构造函数和赋值运算符，防止浅拷贝导致引用失效
    thread_match_tr(const thread_match_tr&) = delete;
    thread_match_tr& operator=(const thread_match_tr&) = delete;
};
extern double is_near_integer_bar ;
void match_TRs(vector<vector<uint32_t>> p, vector<vector<double>> ratios,kmer_TRs_interval& ktis, int kmer_len);
void match_interval(vector<vector<uint32_t>> p, vector<vector<double>> ratios, int kmer_len,vector<pair<int, int>>& interval,vector<int>& interval_space);
void match(char* file_name, int kmer_len, int thread_num, char* distination_file, bool repeated_file_flag , bool overwrite_file=1);
void generate_seq_interval(char* seq, pair<int, int> interval, char** result1);
#endif /* INCLUDE_MATCH_H_ */
