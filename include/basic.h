
#ifndef BASIC_H_
#define BASIC_H_
#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <map>
#include <vector>
#include <algorithm>
#include "stdio.h"
#include <fstream>
#include <ctime>
#include "stdlib.h"
#include <math.h>
#include <iostream>
#include <bitset>
#include <assert.h>
#include <pthread.h>
#include <mutex>
#include <sys/time.h>
#include <unistd.h> 	//access header file
#include <sys/types.h>	//mdkir header file
#include <sys/stat.h> 	//mdkir header file
#include <vector>
#include <cstdint>
#include <cmath>
#include <vector>
#include <fstream>
#include <random>
#include <iomanip>
#include <semaphore.h>// 信号量
#include <utility>
#include <memory>
#include <math.h>
#include <utility> //swap
#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;
using namespace std;

#define radio_bar 0.8

struct bit256KmerPara {
	uint32_t kmer1Len;//k-mer 的总位数
	uint32_t kmer64Len;//64 位块的数量
	uint32_t remainer1to64;//最后一个 64 位块中实际使用的位数
	uint64_t codefor1;
};
struct RefFilePath
{
	char** pRefFilePath;
	uint32_t NumberOfPathes;
};

struct Buffer{
	char* seq=NULL;
	vector<int> pos;
	int length=0;
};
extern Buffer buffer;
extern vector<pair<int,int>>seq_bound;
extern pthread_mutex_t mutex_index;
extern int len_global;
extern int time_global;
extern double edit_distance_radio_global;
extern int init_seq_time;
extern vector<double> ratios_1_in;
extern vector<double> false_charactor_ratio;
extern string file_name_s;



void initBuffer(Buffer& buffer);
void get_para(struct bit256KmerPara* para1, uint32_t kmer_length);
void ReadSeq(char** seq1, uint32_t* seq_length, const char* p_ref);//读取序列
void ReadSeq_ref(char** seq1, uint64_t* seq_length, char* p_ref);
void cal_hash_value_directly_64bit(char* seq, uint64_t &current,uint32_t k_mer_len);
void cal_hash_value_indirectly_64bit(char* seq, uint64_t &current, uint64_t original, uint32_t k_mer_len);
uint32_t cmp256BitKmer(uint64_t* a, uint64_t* b, uint32_t len);
void getRefFilePathes(char* pathFile, struct RefFilePath& p);
uint32_t cmp64BitKmer(uint64_t p, uint64_t b);
void int_value_to_kmer(uint64_t kmer_int, char* kmer, uint32_t kmer_len);
void freeRefFilePath(struct RefFilePath* ref);
void ReadSeq_ref(string& seq, string& filename);

//void init_kti(kmer_TRs_interval& k, int edge_len);


#endif /* BASIC_H_ */
