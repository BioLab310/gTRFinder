
#include "DBGindex.h"
#include "parallel.h"
#include "generateSeq.h"
#include "test.h"
using namespace std;

int main(int argc, char * argv[]) {
	char* example_path=NULL;
	char* init_seq_file=NULL;
	char* distination_file=NULL;
	uint32_t kmer_len=3;
	uint32_t thread_num=3;
	int edge_len=0;
    int len_low = 10;
	int len_up = len_low;
	vector<int> len_lows;
	vector<int> len_ups;
    int times_low = 10;
    int times_up = 10;
    double edit_distance_up = 10;
	double edit_distance_low = 10;
    int seq_num=10;
    int cell_num=1;
    int func_flag=-1;

    //状态标志位置
    bool have_cell_num=0;
    bool have_seq_num=0;
    bool have_init_seq=0;
    bool have_len=0;
    bool have_times=0;
    bool have_ed=0;
    int have_flag_num=6;
	bool repeated_file_flag=false;
    vector<bool> have_flags(6,0);
	for (int i = 0; i < argc; ++i) {
		string arg = argv[i];
		if (arg == "-Sf") {         // graph file
			example_path= argv[++i];
			cout<<argv[i]<<endl;
		}
		else if(arg=="-Of"){
			have_init_seq=1;
			init_seq_file=argv[++i];
			have_flags[0]=1;
		}
		else if(arg== "-kmer_len"){
			kmer_len=atoi(argv[++i]);
			edge_len=kmer_len+1;
			cout<<kmer_len<<endl;
		}
		else if(arg=="-thread_num"){
			thread_num=atoi(argv[++i]);
		}
		else if(arg=="-T"){
			have_times=1;
			times_low=atoi(argv[++i]);
			times_up=atoi(argv[++i]);
			if(times_up ==0 ){
				times_up= times_low;
				i--;
			}
			have_flags[4]=1;
		}
		else if(arg=="-L"){
			have_len=1;
			while(atoi(argv[1+i])!=0){
				len_low=atoi(argv[++i]);
				len_up=atoi(argv[++i]);
				if(len_up ==0 ){
					len_up= len_low;
					i--;
				}
				len_lows.push_back(len_low);
				len_ups.push_back(len_up);
				cout<<len_lows[len_lows.size()-1]<<endl;
			}
			cout<<"len_lows_size"<<len_lows.size()<<endl;
			have_flags[3]=1;
		}
		else if(arg=="-Edr"){
			have_ed=1;
			edit_distance_low=atof(argv[++i]);
			edit_distance_up=atof(argv[++i]);
			if(edit_distance_up ==0.0){
				edit_distance_up= edit_distance_low;
				i--;
			}
			cout<<edit_distance_up<<endl;
			have_flags[5]=1;
		}
		else if(arg=="-Sn"){
			seq_num=atoi(argv[++i]);
			have_seq_num=1;
			have_flags[1]=1;
		}
		else if(arg=="-Cn"){
			cell_num=atoi(argv[++i]);
			have_cell_num=1;
			have_flags[2]=1;
		}
		else if(arg=="-G"){
			func_flag=0;
		}
		else if(arg=="-M"){
			func_flag=1;
		}else if(arg=="-test"){
			func_flag=2;
		}
		else if(arg=="-Df"){
			distination_file = argv[++i];
			// cout<<distination_file<<endl;
		}else if(arg== "-rf"){
			repeated_file_flag= true;
		}else if(arg == "-is_near_integer_bar"){
			is_near_integer_bar=atof(argv[++i]);
			cout<<"is_near_integer_bar"<<' '<<is_near_integer_bar<<endl;
		}
		else {
			// cout<<"error"<<endl;
		}
	}
	py::scoped_interpreter guard{};
	if(func_flag==0){
		generate(len_lows, len_ups, times_low, times_up,cell_num , seq_num, edit_distance_low, edit_distance_up,
		init_seq_file, have_flags, distination_file);
	}else if(func_flag==1){
		//根据文件输入修改
		match(example_path, 7, thread_num , distination_file, repeated_file_flag);
	}else if(func_flag==2){
		test_integer_bar(example_path, 7, thread_num , distination_file, repeated_file_flag);
	}
//
//	struct timeval tvs, tve;
//	double time_token=0;
//
////	gettimeofday(&tvs, NULL);
////	index_gen_seq(seq, seq_len,edge_len);
////	gettimeofday(&tve,NULL);
////	time_token= (tve.tv_sec-tvs.tv_sec)*1e6;
////	time_token= (time_token+ tve.tv_usec- tvs.tv_usec)*1e-6;
////	cout<< " no_parallel take time  "<<time_token<<endl;
//
//	gettimeofday(&tvs, NULL);
//	index_gen_parallel_seq(seq,seq_len, edge_len, thread_num);
//	gettimeofday(&tve,NULL);
//	time_token= (tve.tv_sec-tvs.tv_sec)*1e6;
//	time_token= (time_token+ tve.tv_usec- tvs.tv_usec)*1e-6;
//	cout<< " parallel take time  "<<time_token<<endl;

//	free(seq);
	return 0;
}
