
#include "test.h"
int kmer_radio=3;

void test_integer_bar(char* p_ref, int kmer_len, int thread_num, char* distination_file, bool repeated_file_flag){
	bool overwrite_flag= false;
	for(is_near_integer_bar=0.02; is_near_integer_bar<0.4; is_near_integer_bar+=0.02){
		if(is_near_integer_bar==0.02){
			overwrite_flag = 1;
		}else{
			overwrite_flag = 0;
		}
		match(p_ref, 7, thread_num , distination_file, repeated_file_flag, overwrite_flag);
	}
}
double generateRandom(double lower, double upper) {
    double value= (rand() / (double)RAND_MAX) * (upper - lower) + lower;
	if (value > 100) {
		return (value + 50) / 100 * 100; // 向最近的 100 取整
	} else if (value > 10) {
		return (value + 5) / 10 * 10; // 向最近的 10 取整
	}
//	cout<<value<<endl;
	return value; // 小于或等于 10 的数字保持不变
}

void test_generate(){
	//根据seq_num生成生成lens, times, edit_distance, insert_pos

	srand(time(0)); // 随机种子
	for(int i=1; i<100;i++){
		file_name_s="fasta/"+to_string(i)+".fasta";
		int seq_num=generateRandom(3, 10);
		vector<int> lens;
		vector<int> times(seq_num);
		vector<double> edit_distance;
		vector<int> insert_poss;
		vector<tuple<int, int, int>> results;
		int init_seq_len=5*1e4;
		int insert_pos=0;
		for(int i=0; i< seq_num;i++){
			if(i==0) lens.push_back(generateRandom(10, 100));
			else if(i==1) lens.push_back(generateRandom(2, 6));
			else if(i==2) lens.push_back(generateRandom(100, 300));
			else if(i==3) lens.push_back(generateRandom(300, 1000));
			else lens.push_back(generateRandom(10, 100));

			times[i]=generateRandom(60, 100);
			edit_distance.push_back(generateRandom(0.0,0.10));
			insert_pos=generateRandom(insert_pos,init_seq_len);
			insert_poss.push_back(insert_pos);
		}
		generate_match(init_seq_len,lens, times, edit_distance,insert_poss);
	}
}

void generate_match(int init_seq_len,vector<int> lens, vector<int> times, vector<double> edit_distance, vector<int> insert_pos){
	initBuffer(buffer);
	//生成seq部分
	int kmer_len=7;
	uint32_t seq_len=0;
	char* seq=NULL;
	generate_seq_list_ns(init_seq_len, insert_pos, lens, times, edit_distance, &seq, seq_len);

	//查找匹配部分
	int edge_len=kmer_len+1;
	int thread_num=5;
	vector<vector<uint32_t>> p_int;
	index_gen_parallel_seq1(seq,seq_len, edge_len, thread_num,p_int);
	vector<vector<double>> ratios;
	character_cmp(p_int, ratios);
	cout<<ratios[0].size()<<"识别开始"<<endl;
	vector<pair<int, int>> interval;
	vector<int> interval_space;
	match_interval(p_int, ratios, edge_len,interval,interval_space);
	cout<<"识别结束"<<endl;

	//写入数据
	char* result=NULL;
	const char* file_name=file_name_s.c_str();
	for(int i=0;i< interval.size();i++){
		int cell_length=interval_space[i];
		char write_buffer[100];
		sprintf(write_buffer, "> %d \t%d \t%d \t%d \n", i, cell_length, interval[i].first,interval[i].second );
		printf("> %d \t%d \t%d \t%d  \n", i, cell_length, interval[i].first,interval[i].second);
		result=NULL;
//		sprintf(write_buffer, "> %d \t", i );
		write_file(write_buffer,file_name);
		generate_seq_interval(seq, interval[i], &result);
		write_file(result, file_name);
		char* c="\n";
		write_file(c, file_name);
		result=NULL;
	}
//	free(result);
}

void test_match(int thread_num){
	//生成测试字符串
	char* seq=NULL;
	uint32_t seq_len=0;
	int kmer_len=7;
	vector<int> lens={50,50,50};
	vector<int> times={100,50,50};
	vector<double> edit_distance_radios={0.0,0.05,0.05};
	int init_seq_len=100*lens[0]*times[0];
	vector<int> insert_pos={int(init_seq_len*0.5),int(init_seq_len*0.6),int(init_seq_len*0.7)};
//	vector<int> lens={50};
//	vector<int> times={100};
//	vector<double> edit_distance_radios={0.10};
//	int init_seq_len=100*lens[0]*times[0];
//	vector<int> insert_pos={int(init_seq_len*0.4)};
	generate_seq_list_ns(init_seq_len,insert_pos, lens, times, edit_distance_radios, &seq, seq_len);

	const char *file_name = "lens_50_40_60_times_100_50_20.txt";
	write_file(seq, file_name);

	//区间20000-25000 35000-37000 47000-47600
	//生成p_int ratios
	int edge_len=kmer_len+1;
	vector<vector<uint32_t>> p_int;
	index_gen_parallel_seq1(seq,seq_len, edge_len, thread_num,p_int);
	vector<vector<double>> ratios;
	character_cmp(p_int, ratios);
//	for(int i=0;i<spaces.size();i++){
//		if(spaces[i]>20&& spaces[i]<55){
//			cout<<spaces[i]<<' ';
//		}
//	}
//	kmer_TRs_interval ktis;
	cout<<ratios[0].size()<<"识别开始"<<endl;
//	match_TRs(p_int ,ratios, ktis, edge_len);
	vector<pair<int, int>> interval;
	vector<int> interval_space;
	match_interval(p_int, ratios, edge_len,interval,interval_space);
	cout<<"识别结束"<<endl;
	//识别生成
}

//test_for_edit
void test_histogram( int thread_num){
	int max_kmer_len=9;
	vector<int> kmer_num(10,0);
	int test_seq_num=3;
	vector<int> lens={40,50};
	vector<int> times={50,80,100,120,150,200};
//	vector<double> edit_distance_radios={0.01,0.05,0.1,0.15};
//	vector<int> insert_poss={90000};
//	int init_seq_len=200000;
	for(auto len:lens){//重复单元长度40,50
		len_global=len;
		vector<int> len_test(1, len);
		for(int time=50;time<200;time+=20){//重复单元次数500+200
			time_global=time;
			vector<int> time_test(1, time);
			double max_ed=0.3;
			double add_ed=0.02;
			for(double ed=0.01; ed<max_ed;ed+=add_ed){//编辑距离0.01+0.02
				edit_distance_radio_global=ed;
				vector<double> ed_test(1, ed);
				for(int i=0;i<=30;i=i+3){//总长度1+3
					init_seq_time=i+1;
					int init_seq_len=i*len*time;
					vector<int> ip_test(1, init_seq_len*0.4);
					char* seq = NULL;
					uint32_t seq_len=0;
					generate_seq_list_ns(init_seq_len,ip_test, len_test, time_test, ed_test, &seq, seq_len);
					uint32_t kmer_len=7;//kmer选择
					uint32_t edge_len1=kmer_len+1;
					index_gen_parallel_seq_plot(seq,seq_len, edge_len1, thread_num, 2);
					//开始绘制图形
//						cout<<i+1<<'\t'<<ed<<'\t'<<kmer_len<<endl;
				}
			}
			if(fabs(fmod(max_ed+0.01, add_ed))>1e-6){
				edit_distance_radio_global=max_ed;
				vector<double> ed_test(1, max_ed);
				for(int i=0;i<=30;i=i+3){//总长度1+3
					init_seq_time=i+1;
					int init_seq_len=i*len*time;
					vector<int> ip_test(1, init_seq_len*0.4);
					char* seq = NULL;
					uint32_t seq_len=0;
					generate_seq_list_ns(init_seq_len,ip_test, len_test, time_test, ed_test, &seq, seq_len);
					uint32_t kmer_len=7;//kmer选择
					uint32_t edge_len1=kmer_len+1;
					index_gen_parallel_seq_plot(seq,seq_len, edge_len1, thread_num, 2);
					//开始绘制图形
//						cout<<i+1<<'\t'<<ed<<'\t'<<best_kmer<<endl;
				}
			}
		}
	}

}

void test_kmer( int thread_num){
	int max_kmer_len=9;
	vector<int> kmer_num(10,0);
	int test_seq_num=3;
	vector<int> lens={50};
	vector<int> times={150};
	vector<double> edit_distance_radios={0.01,0.05,0.1,0.15};
//	vector<int> insert_poss={90000};
//	int init_seq_len=200000;
	for(auto len:lens){
		vector<int> len_test(1, len);
		for(auto time:times){
			vector<int> time_test(1, time);
			for(auto ed:edit_distance_radios){
				edit_distance_radio_global=ed;
				vector<double> ed_test(1, ed);
				int best_kmer_of_time=3;
				for(int i=0;i<30;i=i+1){
					init_seq_time=i+1;
					int init_seq_len=i*len*time;
					vector<int> ip_test(1, init_seq_len*0.4);
					char* seq = NULL;
					uint32_t seq_len=0;
//					cout<<"seq_begin"<<endl;
					generate_seq_list_ns(init_seq_len,ip_test, len_test, time_test, ed_test, &seq, seq_len);
					uint32_t kmer_len=3;
					double temp=0;
					int best_kmer=3;
					ratios_1_in.clear();
					false_charactor_ratio.clear();
					for(kmer_len=3;kmer_len<=max_kmer_len;kmer_len++){
						uint32_t edge_len1=kmer_len+1;
						kmer_radio=kmer_len;
						index_gen_parallel_seq(seq,seq_len, edge_len1, thread_num);
					}
					best_kmer=evaluateScore(ratios_1_in, false_charactor_ratio);
//					cout<<best_kmer<<endl;
//					if(best_kmer!=best_kmer_of_time){
					for(kmer_len=3;kmer_len<=max_kmer_len;kmer_len++){
						uint32_t edge_len1=kmer_len+1;
//							cout<<"开始"<<endl;
						kmer_radio=kmer_len;
						index_gen_parallel_seq_plot(seq,seq_len, edge_len1, thread_num, kmer_len==best_kmer);
//							cout<<"结束"<<endl;
//							if(radio_1_in>temp){
//								temp=radio_1_in;
//								best_kmer=kmer_len;
//	//							cout<<temp<<' '<<kmer_len<<endl;
//							}
				//			cout<<len<<' '<<time<<' '<<edit_distance_radio<<' '<< kmer_len<<' '<< radio_1_in<<endl;
					}
					best_kmer_of_time=best_kmer;
					//开始绘制图形
					cout<<i+1<<'\t'<<ed<<'\t'<<best_kmer<<endl;

				}
			}
		}
	}
//	vector<int> lens={50};
//	vector<int> times={1000};
//	vector<double> edit_distance_radios={0.1};
//	vector<int> insert_poss={900000};
//	int init_seq_len=2000000;
//	char* seq = NULL;
//	uint32_t seq_len=0;
//	generate_seq_list_ns(init_seq_len,insert_poss, lens, times, edit_distance_radios, &seq, seq_len);
//	uint32_t kmer_len=3;
//	double temp=0;
//	for(;kmer_len<=max_kmer_len;kmer_len++){
//		uint32_t edge_len1=kmer_len+1;
////		cout<<"开始"<<endl;
//		kmer_radio=kmer_len;
//		index_gen_parallel_seq(seq,seq_len, edge_len1, thread_num);
////		cout<<"结束"<<endl;
//		if(radio_1_in>temp){
//			temp=radio_1_in;
//		}
////			cout<<len<<' '<<time<<' '<<edit_distance_radio<<' '<< kmer_len<<' '<< radio_1_in<<endl;
//	}
//	kmer_num[kmer_radio-1]++;
//	free(seq);
//	 auto maxElementIt = max_element(kmer_num.begin(),kmer_num.end());
//	 int best_kmer_len = distance(kmer_num.begin(), maxElementIt);
//	 cout<<"best_kmer_len"<<best_kmer_len<<endl;
}






