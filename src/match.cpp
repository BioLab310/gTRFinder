
#include "match.h"


bool mind_bar(double n){
	if(n<1){
		n=1/n;
	}
	int x=round(n);
//	cout<<x<<' '<<n<<endl;
	if(abs(x-n)<0.3){
		return true;
	}
	else return false;
}

pthread_mutex_t kti_mutex;
void insert_KTI(kmer_TRs_interval& ktis, char* kmer, startEnd start_end, int kmer_len,int space){
	pthread_mutex_lock(&kti_mutex);
 	// ktis.kmer=(char**)realloc(ktis.kmer, (ktis.length + 1) * sizeof(char*));
 	// ktis.kmer[ktis.length]=(char*)malloc(sizeof(char)*kmer_len+1);
 	// strcpy(ktis.kmer[ktis.length],kmer);
 	ktis.length++;
 	ktis.start_ends.push_back(start_end);
 	ktis.spaces.push_back(space);
	pthread_mutex_unlock(&kti_mutex);

	// if(ktis.length){
	// 	cout<<kmer<<'\t';
	// 	for(int i=0;i<start_end.size();i++){
	// 		cout<<start_end[i].first<<' '<<start_end[i].second<<'\t';
	// 	}
	// 	cout<< space<<endl;
	// }

}
int compute_space(vector<uint32_t> p){
	double epsilon=0.1;
	int threshold=5;
	if(p.size()<=1) return 1e5;
	vector<int> pds(p.size()-1);
//	cout<<pds.size()<<endl;
	for(int i=0;i< pds.size();i++){
		pds[i]=p[i+1]-p[i];
	}
	int32_t min_ratio=1e7;
	int space=*max_element(pds.begin(), pds.end());
	for(int i=0;i< pds.size();i++){
		int j=i+1;
		int temp_space=pds[i];
		int count=1;
		for(;j<pds.size();j++){
			if(abs(pds[j]-temp_space)< temp_space*epsilon){
				count++;
			}
		}
		if(count>=threshold and temp_space<min_ratio){
			space=temp_space;
			min_ratio=temp_space;
		}
	}
//	cout<<"差值列表"<<endl;
//	for(auto i:pds){
//		cout<<i<<' ';
//	}
//	cout<<endl;
//	cout<<space<<endl;
	return space;

}
void space_push_back(vector<int>& spaces, int space, int len){
	for(int i= 0; i< len; i++){
		spaces.push_back((space));
	}
}
void sort_interval_space(vector<pair<int, int>>& interval, vector<int>& interval_space){
	// 创建索引数组
    vector<int> indices(interval.size());
    for (size_t i = 0; i < indices.size(); i++) {
        indices[i] = i;
    }

    // 对索引数组进行排序，并通过索引同步排序 interval 和 space
    sort(indices.begin(), indices.end(), [&](int i, int j) {
        return interval[i].first < interval[j].first;
    });

    // 使用排序后的索引重构 interval 和 space
    vector<pair<int, int>> sorted_interval;
    vector<int> sorted_space;
    for (int idx : indices) {
        sorted_interval.push_back(interval[idx]);
        sorted_space.push_back(interval_space[idx]);
    }

    // 替换原有的数据
    interval = move(sorted_interval);
    interval_space = move(sorted_space);
}
void add_score(double &score, int add_score){
	if(add_score> 1.5){
		score+= 1.0/add_score;
	}else if(add_score >0.75){
		score+=add_score;
	}else {
		score += add_score;
	}
}
double is_near_integer_bar = 0.3;
bool isNearInteger(double value) {
	if(value < 1) {
		value =1/value;//取倒数
	}
    // 获取最接近的整数
    double nearestInteger = std::round(value);
    // 检查 value 与 nearestInteger 之间的差值
    return std::abs(value - nearestInteger) <= is_near_integer_bar;
}




void match_TRs(vector<vector<uint32_t>> p, vector<vector<double>> ratios,kmer_TRs_interval& ktis, int kmer_len ){

	int ii=0;
	double bar= 6;
	for(int i=0;i<ratios.size();i++){
		if(ratios[i].size()<3) continue;
		char* kmer=(char*)malloc(sizeof(char)*kmer_len);
		int_value_to_kmer(i, kmer, kmer_len);
		startEnd start_end;
		int start=0;
		int end=0;
		int space= compute_space(p[i]);
		double score=0;
		double score_bar= 2;
		int continue_time = 0;
		for(int j=0;j<ratios[i].size();j++){
			// cout<< kmer<< '\t'<< j<< '\t'<<ratios.size()<<endl;
			//增加一个条件 向后扩展的区间范围不能大于5倍的space
			add_score(score, ratios[i][j]);
			bool continue_flag= ratios[i][j]<bar&&ratios[i][j]>1/bar&&(p[i][end]-p[i][end-1])<=10*space&& isNearInteger(ratios[i][j]);
			if(continue_flag){
				end=j+1;
			}
			else{
				int TR_len=end-start;
				if(score>score_bar){
					start_end.emplace_back(p[i][start], p[i][end]);
				}
				start=j+1;
				end=start;
				score = 0 ;
			}
		}

		int TR_len=end-start;
		if(score>score_bar){
			start_end.emplace_back(p[i][start], p[i][end]);
		}
		if(start_end.size()>0){

			insert_KTI(ktis, kmer ,start_end,kmer_len, space);
		}
		else{
		}
		free(kmer);
	}
}

vector<pair<int,int>> merge_interval(vector<pair<int,int>> interval_src, vector<int>& interval_space_result){
	vector<pair<int, int>> interval_result;
	if(interval_src.size()<1) return interval_result;
	vector<int> interval_space = interval_space_result;
	interval_space_result.clear();
	sort_interval_space(interval_src, interval_space);
	
	//[a, b] [c,d]区间有重叠部分,组合为[a,d]
	int start=interval_src[0].first;
	int end=interval_src[0].second;
	int i=0;
	int interval_sum=0;
	int space_sum=0;
	for(i=0;i<interval_src.size()-1;i++){
		if(end>=interval_src[i+1].first){
			start=start;
			end=max(interval_src[i+1].second, end);
			interval_sum++;
			space_sum+= interval_space[i];
		}
		else{
			if(interval_sum>0){
				interval_result.emplace_back(start, end);
				interval_space_result.push_back(space_sum/interval_sum);
				// cout<<space_sum<<' '<< interval_sum<<' '<<space_sum/interval_sum<<endl;
			}
			space_sum=0;
			interval_sum=0;
			start=interval_src[i+1].first;
			end=interval_src[i+1].second;
		}
	}
	if(interval_sum>0){
		interval_result.emplace_back(start, end);
		interval_space_result.push_back(space_sum/interval_sum);
		// cout<<space_sum<<' '<< interval_sum<<' '<<space_sum/interval_sum<<endl;
	}
	return interval_result;
}

// 线程执行的函数
void* thread_match_TRs(void* arg) {
    thread_match_tr* args = (thread_match_tr*)arg;

    // 在这里处理分配给该线程的子集
    vector<vector<uint32_t>>& p = args->p_subset;
    vector<vector<double>>& ratios = args->ratios_subset;
    kmer_TRs_interval& ktis = args->ktis;
    int kmer_len = args->kmer_len;
	match_TRs(p ,ratios, ktis, kmer_len);//生成kmer与对应的区间
	delete args;
    pthread_exit(NULL); // 线程退出
}

void merge_ktls(kmer_TRs_interval ktis, vector<pair<int,int>>& result, int kmer_length, vector<int>& interval_space){
	vector<pair<int, int>> temp_interval;
	vector<int> temp_interval_space;
	ktis.sortBySpace();
	vector<pair<int, int>> interval;
	if(ktis.length<=0) return ;
//	cout<<"length"<<ktis.length<<endl;
	int temp_space=ktis.spaces[0];
	int first_space=temp_space;
	int kmer_sum=0;
	for(int i=0;i<ktis.length;i++){
		// space与当前kmer的space大致相同,加入合并区间
		bool is_over_space = abs(temp_space - first_space) > min(2,int(first_space*0.01));
		if(ktis.spaces[i]>temp_space*1.1 || is_over_space){
			// cout<<first_space<<' '<<temp_space<<' '<<kmer_sum<<endl;
			if(kmer_sum>=temp_space*0.1){
//				for(auto j:temp_interval){
//					cout<<j.first<<' '<<j.second<<endl;
//				}
				interval=merge_interval(temp_interval, temp_interval_space);
				result.insert(result.end(), interval.begin(), interval.end());//加入结果中
				interval_space.insert(interval_space.end(), temp_interval_space.begin(), temp_interval_space.end());
			}
			temp_space=ktis.spaces[i];//初始化space大小相近的区间
			first_space=ktis.spaces[i];
			temp_interval.clear();
			temp_interval_space.clear();
			kmer_sum=0;
		}
		else{
			temp_interval.insert(temp_interval.end(),ktis.start_ends[i].begin(),ktis.start_ends[i].end());
			space_push_back(temp_interval_space, ktis.spaces[i], ktis.start_ends[i].size());
//			cout<<ktis.start_ends[i].size()<<' '<<ktis.start_ends[i][0].first<<endl;
			temp_space=ktis.spaces[i];
			kmer_sum++;
		}
	}
	if(kmer_sum>=temp_space*0.1){
		interval=merge_interval(temp_interval, temp_interval_space);
		result.insert(result.end(), interval.begin(), interval.end());//加入结果中
		interval_space.insert(interval_space.end(), temp_interval_space.begin(), temp_interval_space.end());
	}
	// for(int j=0;j<result.size();j++){
	// 	cout<<j+1<<'\t'<<result[j].first<<'-'<<result[j].second<<'-'<<interval_space[j]<<' ';
	// }
	// cout<<endl;
}


void overlap_inverse_merge(vector<pair<int, int>>& interval, vector<int>& interval_space){
	if(interval.size()<2) return;
	
	sort_interval_space(interval,interval_space);

	// 检查[a,b] [c, d] [e,f]之间是否有重叠
	for (int i = 0; i < interval.size() - 1; ) {
	    int j = i + 1;
//	    cout<< interval_space[i] <<'<'<< interval_space[j]<<endl;
	    while (j < interval.size() && interval[i].second> interval[j].first+10) {
	        bool space_j_more = interval_space[i] > interval_space[j];
	        if (!space_j_more) {
	            interval.erase(interval.begin() + j);
	            interval_space.erase(interval_space.begin() + j);

	        } else {
	            interval.erase(interval.begin() + i);
	            interval_space.erase(interval_space.begin() + i);

	            break;
	        }

	    }
	    // 如果没有删除 `i`，则增加 `i`
//	    cout<<i<<' '<<j<<endl;
	    if (j >= interval.size() || interval[i].second <= interval[j].first+10) {
	        i++;
	    }
	}
}


void match_interval(vector<vector<uint32_t>> p, vector<vector<double>> ratios, int kmer_len,vector<pair<int, int>>& interval,vector<int>& interval_space){
	kmer_TRs_interval ktis;
	ktis.length=0;
	ktis.kmer=NULL;
	int p_size = p.size();
    int ratios_size = ratios.size();
	int num_threads = 100;
	int chunk_size = p_size / num_threads;
    int remainder = p_size % num_threads;
	// 2. 创建线程
    // pthread_t threads[num_threads];
	// pthread_mutex_init(&kti_mutex, NULL);
    // // 3. 分配数据给每个线程
    // int start_index = 0;
    // for (int i = 0; i < num_threads; ++i) {
    //     int current_chunk_size = min( chunk_size + (i+1 < remainder ? 1 : 0),int( p.size())); // 均匀分配余数
	// 	thread_match_tr* arg = new thread_match_tr(ktis,
    //     vector<vector<uint32_t>>(p.begin() + start_index, p.begin() + start_index + current_chunk_size),
    //     vector<vector<double>>(ratios.begin() + start_index, ratios.begin() + start_index + current_chunk_size),
    //     kmer_len); 
    //     // 创建线程
    //     int rc = pthread_create(&threads[i], NULL, thread_match_TRs, (void*)arg);
    //     if (rc) {
    //         cerr << "Error: pthread_create() failed, return code is " << rc << endl;
    //         exit(-1);
    //     }

    //     start_index += current_chunk_size;
    // }
	// for (int i = 0; i < num_threads; ++i) {
    //     pthread_join(threads[i], NULL);
    // }
	match_TRs(p ,ratios, ktis, kmer_len);
	interval_space.resize(0);//区间对应的space
	merge_ktls(ktis, interval, kmer_len, interval_space);//对区间进行合并
//	cout<<"merge_TRs"<<endl;
	overlap_inverse_merge(interval, interval_space);
	pthread_mutex_destroy(&kti_mutex);
	// for(int j=0;j<interval_space.size();j++){
	// 	cout<<j+1<<'\t'<<interval_space[j]<<' ';
	// }
//	cout<<endl;

}

void generate_seq_interval(char* seq, pair<int, int> interval, char** result1){
	char* result;
	if(interval.first> interval.second||interval.first<0) return ;
	int len=interval.second-interval.first+1;
	result=(char*) malloc (sizeof(char)*(len+1));
	if(interval.second>strlen(seq)) return ;
	strncpy(result, seq+interval.first, len );
	result[len]='\0';
	*result1=result;
}

string extract_filename(char* file_name) {
	const char* last_slash = strrchr(file_name, '/');
	const char* base_name = (last_slash != nullptr) ? (last_slash + 1) : file_name;

	// 找到最后一个点
	const char* last_dot = strrchr(base_name, '.');
	if (last_dot != nullptr) {
		// 去掉后缀名
		return std::string(base_name, last_dot - base_name);
	}
	// 如果没有点，返回整个文件名
	return std::string(base_name);
}
string generate_match_filename(const string& base_path, bool repeated_file_flag) {
    // 构造初始文件路径
    std::string file_name = base_path+".fasta" ;
    std::ifstream file(file_name);
    int counter = 1;
	if(!repeated_file_flag) return file_name;
    // 如果文件已存在，生成新的路径
    while (file.good()) {
        file_name = base_path+"(" + std::to_string(counter) +")"+ ".fasta";
        file.close();
        file.open(file_name);
        ++counter;
    }
    return file_name;
}

struct match_data{
	int kmer_len;
	int thread_num;
	char* seq;
	uint64_t seq_len;
	int seq_count;
	char* match_file_name;
	// 构造函数
    match_data(int kmer_len, int thread_num, const char* seq_, uint64_t seq_len, int seq_count, const char* match_file_name_)
        : kmer_len(kmer_len),
          thread_num(thread_num),
          seq_len(seq_len),
          seq_count(seq_count)
    {
		seq = new char[seq_len+1];
		strncpy(seq, seq_, seq_len);
		match_file_name = new char[strlen(match_file_name_) + 1];
		strcpy(match_file_name, match_file_name_);
    }

    // 析构函数
    ~match_data() {
        delete[] seq;             // 释放 seq 的内存
        delete[] match_file_name; // 释放 match_file_name 的内存
		seq = NULL; 
		match_file_name= NULL;
    }
};

pthread_mutex_t file_mutex;
pthread_mutex_t seq_mutex;
sem_t match_semaphore;
//  写入文件锁的初始化
void match_write(int kmer_len, int thread_num, char *seq, uint64_t seq_len, int seq_count,const char* match_file_name){
	int edge_len=kmer_len+1;
	//生成序列对应的位置列表
	vector<vector<uint32_t>> p_int;
	// 创建多线程识别
	thread_num = min(int(seq_len / 100000) +1 , 100);
	index_gen_parallel_seq1(seq,seq_len, edge_len, thread_num,p_int);
	//分析匹配数据
	vector<vector<double>> ratios;
	character_cmp(p_int, ratios);
	vector<pair<int, int>> interval;
	vector<int> interval_space;
	match_interval(p_int, ratios, edge_len,interval,interval_space);

	//将结果写入fasta文件中
	char* result=NULL;
	for(int i=0;i< interval.size();i++){
		int cell_length=interval_space[i];
		pthread_mutex_lock(&file_mutex);
		char write_buffer[100];
		sprintf(write_buffer, "> %d \t%d \t%d \t%d \t%d \t%.2f\n",seq_count+1, i+1, cell_length, interval[i].first,interval[i].second,is_near_integer_bar);
		// printf("> %d \t%d \t%d \t%d \t%d \t%f\n",seq_count+1, i+1, cell_length, interval[i].first,interval[i].second,is_near_integer_bar);
		result=NULL;
//		sprintf(write_buffer, "> %d \t", i );
		write_file(write_buffer,match_file_name);
		pthread_mutex_unlock(&file_mutex); 

		generate_seq_interval(seq, interval[i], &result);
		pthread_mutex_lock(&file_mutex);

		write_file(result, match_file_name);
		char* c="\n";
		write_file(c, match_file_name);
		pthread_mutex_unlock(&file_mutex);
		
		free(result);
		result=NULL;
	}

}
void* thread_match(void *arg){
	match_data* a= (match_data*)arg;
	match_write(a->kmer_len, a-> thread_num, a->seq, a->seq_len, a->seq_count, a->match_file_name);
	delete a;
	sem_post(&match_semaphore); //数据释放后释放
	pthread_exit(NULL);
}
int count_to_next_symbol_optimized(FILE* fp, char symbol = '>') {
    const size_t buffer_size = 4096; // 4KB 缓冲区
    char buffer[buffer_size];
    long current_pos = ftell(fp);
    int total_count = 0;

    while (fgets(buffer, buffer_size, fp) != NULL) {
        char* found = strchr(buffer, symbol);
        if (found) {
            // 找到符号，计算位置并返回
            fseek(fp, current_pos, SEEK_SET);
            return total_count + (found - buffer);
        }
        total_count += strlen(buffer);
//        cout<<total_count<<endl;
    }
    // 未找到符号
    fseek(fp, current_pos, SEEK_SET);
    return total_count;
}
void match(char* p_ref, int kmer_len, int thread_num, char* distination_file, bool repeated_file_flag, bool overwrite_file){

	// 打开读出文件，并初始化缓冲区
	uint32_t buffer_size = 1024;
	char buffer_line[1024];
	memset(buffer_line, 0, buffer_size);
	FILE* fp;
	fp = fopen(p_ref,"r+");
	if (fp == NULL)
	{
		cout << "file can not be open!" << endl;
		return;
	}

	// 产生写入文件，返回一个文件打开符fp
	string base_filename = "match_fasta/";
	string file_name_s;
	if(distination_file==NULL) {
		base_filename=base_filename+extract_filename(p_ref);
		file_name_s=generate_match_filename(base_filename, repeated_file_flag);
	}else {
		base_filename= base_filename+ string(distination_file);
		file_name_s= base_filename;
		cout<<"match result file name is: "<<file_name_s<<endl;	
	}
	

	//初始化序列
	char ch;
	char* seq=NULL;
	int total_size=count_to_next_symbol_optimized(fp);
	int last_seq_size=0;
	seq=(char*) malloc (sizeof(char) *total_size*2);
	uint64_t seq_len = 0;
	int seq_count=0;
	const char* match_file_name=file_name_s.c_str();
	if(overwrite_file){
		write_file("", match_file_name,"w");
	}
	int in_sequence = 0; // 0: not in sequence, 1: in sequence
	pthread_mutex_init(&file_mutex, NULL);
	pthread_mutex_init(&seq_mutex, NULL);
	int max_thread= 10;
	sem_init(&match_semaphore, 0, max_thread);
	vector<pthread_t> threads;
	while (fgets(buffer_line, buffer_size - 1, fp) != NULL)
	{
		for (uint32_t i = 0; i < strlen(buffer_line); i++) {
            if (buffer_line[i] == '>') {
                if (in_sequence) {
                    // 结束序列
					sem_wait(&match_semaphore); // 创建数据前加锁
					match_data* threadData = new match_data(kmer_len, thread_num, seq, seq_len, seq_count,  match_file_name);
					pthread_t threadId;
					int rc = pthread_create(&threadId, NULL, thread_match, (void *)threadData);
					threads.push_back(threadId); // 存储线程ID
					// match_write(kmer_len,  thread_num, seq, seq_len, seq_count, match_file_name);
					// cout<<seq_count<<"\t"<<seq_len<<endl;
                    seq_len = 0;
                    seq_count++;
                    in_sequence = 1; // 退出序列读取状态
					// 开始新的序列, 或仅是 '>' 行的开始标记
                    total_size = count_to_next_symbol_optimized(fp); // 计算到下一个 '>' 的大小
                    if (total_size > last_seq_size) {
                        if (seq != NULL) {
                             free(seq); // 释放旧的 seq 缓冲区
                        }
                        seq = (char*)realloc(NULL, sizeof(char) * total_size * 2); // 注意这里使用 NULL 作为 realloc 的第一个参数以进行初始分配
                        if (seq == NULL) {
                            perror("Failed to reallocate seq buffer");
                            fclose(fp);
                            return ;
                        }
                        last_seq_size = 2 * total_size;
                    }else if (seq == NULL) {
                         seq = (char*)malloc(sizeof(char) * last_seq_size); // 如果之前没有分配过 seq，则初始分配
                         if (seq == NULL) {
                            perror("Failed to allocate seq buffer");
                            fclose(fp);
                            return ;
                        }
                    }
                    memset(seq, 0, last_seq_size);  // 初始化 seq 缓冲区
                    seq_len = 0;
                    continue; // 处理完 '>' 后继续处理剩余的字符，或者跳出 inner for 循环处理下一行
                } else {
                    // 开始新的序列, 或仅是 '>' 行的开始标记
                    total_size = count_to_next_symbol_optimized(fp); // 计算到下一个 '>' 的大小
                    if (total_size > last_seq_size) {
                        if (seq != NULL) {
                             free(seq); // 释放旧的 seq 缓冲区
                        }
                        seq = (char*)realloc(NULL, sizeof(char) * total_size * 2); // 注意这里使用 NULL 作为 realloc 的第一个参数以进行初始分配
                        if (seq == NULL) {
                            perror("Failed to reallocate seq buffer");
                            fclose(fp);
                            return ;
                        }
                        last_seq_size = 2 * total_size;
                    } else if (seq == NULL) {
                         seq = (char*)malloc(sizeof(char) * last_seq_size); // 如果之前没有分配过 seq，则初始分配
                         if (seq == NULL) {
                            perror("Failed to allocate seq buffer");
                            fclose(fp);
                            return ;
                        }
                    }
                    memset(seq, 0, last_seq_size);  // 初始化 seq 缓冲区
                    seq_len = 0;
                    in_sequence = 1; // 进入序列读取状态
                    continue; // 处理完 '>' 后继续处理剩余的字符，或者跳出 inner for 循环处理下一行
                }
            } else if (in_sequence) {
                // 在序列中，处理字符
                if (buffer_line[i] >= 'a' && buffer_line[i] <= 'z') {
                    buffer_line[i] -= 32; // 转换为大写
                }
                if (buffer_line[i] != 'A' && buffer_line[i] != 'C' && buffer_line[i] != 'G' && buffer_line[i] != 'T' && buffer_line[i] != '\n' && buffer_line[i] != '\r') { // 保留换行符的判断，如果你的数据可能包含，并需要忽略, 并且移除 '\r'
                    // buffer_line[i] = ''; // 替换为 'A'
					continue;
                }
                if (buffer_line[i] == '\n' || buffer_line[i] == '\r') {
                    continue; // 忽略换行符和回车符
                }
				
                seq[seq_len] = buffer_line[i];
                seq_len++;
                if (seq_len >= total_size*2) { // 理论上不应该超出，但作为安全检查
                    printf("Warning: Sequence length exceeds expected total_size.\n");
                    // break; // 停止当前行的处理，可能需要更精细的错误处理
                }
            }
            // 如果 !in_sequence 且不是 '>'，则忽略或根据需求处理其他非序列数据
        }
		
	}
	if(seq_len>0&& seq!=NULL){
		sem_wait(&match_semaphore); // 创建数据前加锁
		match_data* threadData = new match_data(kmer_len, thread_num, seq, seq_len, seq_count,  match_file_name);
		pthread_t threadId;
		int rc = pthread_create(&threadId, NULL, thread_match, (void *)threadData);
		threads.push_back(threadId); // 存储线程ID
		memset(seq,0,seq_len+1);
		seq_len=0;
		seq_count++;
	}
	for (pthread_t threadId : threads) {
        pthread_join(threadId, NULL); // 等待每个线程完成
    }
	
	sem_destroy(&match_semaphore);// 销毁信号量
	//   销毁互斥锁
    if (pthread_mutex_destroy(&file_mutex) != 0) {
        std::cerr << "Mutex destruction failed" << std::endl;
        return ;
    }
	cout<<"序列数量 "<<seq_count<<endl;
	free(seq); // 确保释放
}
