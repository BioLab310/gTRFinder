

#ifndef INCLUDE_TEST_H_
#define INCLUDE_TEST_H_
#include "basic.h"
#include "parallel.h"
#include "evaluate.h"
#include "generateSeq.h"
#include "match.h"
struct kmer_ed{
	uint32_t kmer;//kmer长度
	uint32_t ed;//编辑距离
};
//extern kmer_ed test_num;
extern int kmer_radio;
void test_histogram(int thread_num);
void test_kmer( int thread_num);
void test_match( int thread_num);
void test_integer_bar(char* p_ref, int kmer_len, int thread_num, char* distination_file, bool repeated_file_flag);

void test_generate();
void generate_match(int init_seq_len,vector<int> lens, vector<int> times, vector<double> edit_distance, vector<int> insert_pos);
#endif /* INCLUDE_TEST_H_ */
