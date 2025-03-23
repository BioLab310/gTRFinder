#ifndef DBGindex_H_
#define DBGindex_H_
#include "basic.h"
#include "hash.h"
#include "generateSeq.h"

struct index{
	int** p;
	int* len;
	int* maxlen;
};
void index_gen_disp(char* pathFile, uint32_t kmer_len);
void generate_index(bit64Hash* p_root, char* seq, uint64_t seq_length, uint32_t kmer_len);
void index_display(bit64Hash* ph, uint32_t kmer_len);
void generate_index_int(vector<vector<uint32_t>>& p, char* seq,uint64_t begin, uint64_t end, uint32_t kmer_len);
void index_gen_seq(char* seq, uint32_t seq_length,uint32_t kmer_len);
void index_display_int(vector<vector<uint32_t>> p, uint32_t kmer_len ,vector<vector<double>> character, int seq_length, int flag);
void character_cmp(vector<vector<uint32_t>> index,vector<vector<double>>& charactor);
#endif

