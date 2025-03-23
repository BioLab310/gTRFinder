

#ifndef GENERATESEQ_H_
#define GENERATESEQ_H_
#include "basic.h"
#include "test.h"
void generate_seq_list(uint32_t start_len, int len, int times, double edit_distance, char** seq11, uint32_t& seq_len) ;
vector<int> find_top_radio_index(vector<vector<double>> data, int seq_length);
void generate_seq_list_edit(int len, int time,double edit_distance_radio, char**seq11,uint32_t& seq_length);
void generate_seq(int length, char** seq1);
void generate_seq_list_ns(int init_seq_len,vector<int> insert_poss, vector<int> lens, vector<int> times,
		vector<double> edit_distance_radios, char** seq11, uint32_t& seq_len,  bool have_init_seq=0, char* init_seq=NULL);
void generate(vector<int> len_lows, vector<int> len_ups, int time_low, int time_up,  int cell_num, int seq_num, 
	double ed_low,double ed_up, char* init_seq_file,vector< bool > have_flags, char* distination_file);
int write_file(char* data, const char* file_name, const char * flag);
int write_file(char* data, const char* file_name);
int generateRandom(int lower, int upper);
#endif /* GENERATESEQ_H_ */
