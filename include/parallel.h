

#ifndef PARALLEL_H_
#define PARALLEL_H_
#include "DBGindex.h"
#include "plot_data.h"
struct thread_data{
	FILE* fp;
	vector<vector<uint32_t>> index;
	char* seq;
	uint32_t begin;
	uint32_t end;
	uint32_t kmer_len;
	thread_data(FILE* rp, vector<vector<uint32_t>> in, char* s, uint32_t bg, uint32_t ed, uint32_t kl):
		fp(rp), index(in), seq(s), begin(bg), end(ed), kmer_len(kl){}
    ~thread_data() {
    }
    thread_data() : fp(nullptr), seq(nullptr), begin(0), end(0), kmer_len(0) {}
};
void* seq_index_parallel(void* arg);
void index_gen_parallel(char* pathFile, uint32_t kmer_len,int thread);
void ReadSeq_parallel(FILE* fp,char* seq, uint64_t begin, uint64_t end);
void index_gen_parallel_seq1(char* seq, uint32_t seq_length, uint32_t kmer_len,int thread, vector<vector<uint32_t>>& p_int);//耦合度降低
void index_gen_parallel_seq(char* seq, uint32_t seq_length, uint32_t kmer_len,int thread );
void index_gen_parallel_seq_plot(char* seq, uint32_t seq_length, uint32_t kmer_len,int thread,int flag);//绘制图形时候使用
#endif /* PARALLEL_H_ */
