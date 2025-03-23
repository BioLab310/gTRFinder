
#include "generateSeq.h"
vector<string> seq_cell;
bool is_repeating_prefix(const string& seq, int prefix_len) {
    if (seq.length() < 2 * prefix_len) {
        return false; // 序列长度不足以构成重复
    }
    return seq.substr(seq.length() - prefix_len) == seq.substr(seq.length() - 2 * prefix_len, prefix_len);
}

void generate_seq(int length, char** seq1) {
    char acgt[4] = { 'A', 'C', 'G', 'T' };
    string seq = ""; // 使用 string 便于操作
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, 3);
	int check_threshold =10;
    for (int i = 0; i < length; ++i) {
        bool valid = false;
        char new_char;
        while (!valid) {
            new_char = acgt[dis(gen)]; // 生成新字符
            valid = true;
            seq += new_char; // 先尝试添加到序列中
			if(i+1 < check_threshold)
            {	// 检查是否存在重复前缀 (长度从 1 到 i/2)
				for (int prefix_len = 1; prefix_len <= (i + 1) / 2; ++prefix_len) {
					if (is_repeating_prefix(seq, prefix_len)) {
						valid = false; // 发现重复，需要重新生成字符
						seq.pop_back(); // 移除刚才添加的字符
						break;
					}
				}
			}
        }
    }

    // 将 string 转换为 char*
    *seq1 = (char*)malloc(sizeof(char) * (length + 1));
    if (!*seq1) {
        cerr << "Memory allocation failed!" << endl;
        exit(EXIT_FAILURE);
    }
    strcpy(*seq1, seq.c_str());
}

void generate_perfect_seq(int len, int times, char** seq1) {
    char* cell = NULL;
    generate_seq(len, &cell);
//    cout<<cell<<endl;
    int length = len * times;
    char* p_list = (char*)malloc(sizeof(char) * (length + 1));
    if (!p_list) {
        cerr << "Memory allocation failed!" << endl;
        free(cell); // 确保在错误情况下释放cell
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < length;) {
        for (int j = 0; j < len; j++) {
            p_list[i] = cell[j];
            i++;
        }
    }
    p_list[length] = 0;
    *seq1 = p_list;
	seq_cell.push_back(string(cell));
	// free(cell); // 释放cell，避免内存泄漏
}

void shuffle(int* array, int n) {
    for (int i = n - 1; i > 0; i--) {
        int j = rand() % (i + 1); // 生成0到i之间的随机数
        // 交换 array[i] 和 array[j]
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
}

void generate_seq_list_ns(int init_seq_len,vector<int> insert_poss, vector<int> lens, vector<int> times,
		vector<double> edit_distance_radios, char** seq11, uint32_t& seq_len, bool have_init_seq, char* init_seq) {
	seq_bound.resize(0);
	if(insert_poss.size() != times.size() ||
		    insert_poss.size() != lens.size() ||
		    insert_poss.size() != edit_distance_radios.size()) {
		cout<<insert_poss.size()<<' '<<times.size()<<' '<<times.size()<<' '<<edit_distance_radios.size()<<endl;
		cout<<"have no equal characters"<<endl;
		return ;
	}
	if(insert_poss[insert_poss.size()-1]>init_seq_len) {
		cout<<"insert_poss>init_seq_len"<<endl;
		return;
	}
	const char *file_name = file_name_s.c_str();

	char* seq=(char*)malloc(sizeof(char)*(1+init_seq_len));
	memset(seq, 0 ,init_seq_len+1);

	//判断是否随机生成初始序列，或者复已有序列
	if(have_init_seq==0){
		generate_seq(init_seq_len,&seq);
	}else {
		strncpy(seq, init_seq, init_seq_len);
	}
	//插入重复序列
	seq_len=init_seq_len+1;
	for(int i = 0; i < lens.size(); i++) {
	    int insert_pos = seq_len - init_seq_len + insert_poss[i]; // 设置插入位置
	    char* seq_temp = nullptr;
	    uint32_t seq_temp_len = 0;
	    generate_seq_list_edit(lens[i], times[i], edit_distance_radios[i], &seq_temp, seq_temp_len);
//	    sprintf(buffer, "> insert_pos :%d-%d \n",
//	    		insert_pos, insert_pos+ seq_temp_len);
//		write_file(buffer, file_name);
//		printf(" > %d %d, %.2f, %d-%d\n",
//	    		lens[i], times[i], edit_distance_radios[i], insert_pos, insert_pos+ seq_temp_len);
	    seq_len += seq_temp_len;
	    // 重新分配内存，确保成功分配
	    char* temp_seq = (char*)realloc(seq, sizeof(char) * seq_len);
	    if (!temp_seq) {
	        free(seq);
	        std::cerr << "Memory allocation failed!" << std::endl;
	        return; // 或其他错误处理
	    }
	    seq = temp_seq;
	    // 从后向前移动元素，避免覆盖
	    for (int j = seq_len - 1; j >= insert_pos + seq_temp_len; --j) {
	        seq[j] = seq[j - seq_temp_len];
	    }
	    // 将 seq_temp 中的元素复制到 seq 中
	    for (int j = 0; j < seq_temp_len; j++) {
	        seq[insert_pos + j] = seq_temp[j];
	    }
	    seq_bound.push_back(make_pair(insert_pos, insert_pos+seq_temp_len));
	    free(seq_temp);
		
	    seq_temp_len = 0;
//	    cout<<seq<<endl;
	}
	*seq11=seq;
//	char* c="\n";
//	write_file(c, file_name);
//	write_file(seq, file_name);
}
void generate_seq_list_edit(int len, int timei,double edit_distance_radio, char**seq11,uint32_t& seq_length){
	uint32_t seq_length1 = timei* len;
	uint8_t* label = (uint8_t*)malloc(sizeof(uint8_t) * seq_length1);
	if (!label) {
		cerr << "Memory allocation failed!" << endl;
		exit(EXIT_FAILURE);
	}
	memset(label, 0, seq_length1 * sizeof(uint8_t));
	random_device rd;
	mt19937 gen(rd());
	int edit_distance = seq_length1 * edit_distance_radio;
	if (seq_length1 < edit_distance) {
		free(label); // 释放内存
		return;
	}

	uniform_int_distribution<> dis_label(1, 3);
	uniform_int_distribution<> dis_insert(0, seq_length1);
	int32_t edit_length = 0;
	int* index = (int*)malloc(sizeof(int) * seq_length1);
	if (!index) {
		cerr << "Memory allocation failed!" << endl;
		free(label);
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < seq_length1; i++) {
		index[i] = i;
	}
	srand(time(NULL)); // 设置随机种子
	shuffle(index, seq_length1);
    for (size_t i = 0; i < edit_distance; i++) {
        label[index[i]] = dis_label(gen); // 随机修改为1, 2, 或3, 1删除，2替换，3插入
        edit_length = edit_length + (label[index[i]] - 2);
    }
	seq_length = edit_length + seq_length1;
	int insert_pos = dis_insert(gen);
	char* seq = (char*)malloc(sizeof(char) * (seq_length + 1));
	if (!seq) {
		cerr << "Memory allocation failed!" << endl;
		free(label);
		free(index);
		exit(EXIT_FAILURE);
	}
	memset(seq, '\0', seq_length + 1);
	char* seq1 = NULL;
	char* temp=NULL;
	generate_perfect_seq(len, timei, &seq1);
//	cout<<"generate_perfect_seq finished"<<endl;
	uint32_t i = 0, j = 0;
	for (; i < seq_length1; i++) {
		if (label[i] == 0) {
			seq[j] = seq1[i];
		} else if (label[i] == 1) {
			continue; // 删除操作
		} else if (label[i] == 2) {
			generate_seq(1, &temp);
			seq[j] = temp[0];
			free(temp); // 释放temp
		} else if (label[i] == 3) {
			generate_seq(1, &temp);
			seq[j] = temp[0];
			j++;
			seq[j] = seq1[i];
			free(temp); // 释放temp
		}
		j++;
	}
	*seq11=seq;
//	free(temp);
	free(seq1);
	free(label);
	free(index);
}
void generate_seq_list(uint32_t start_len, int len, int times, double edit_distance_radio, char** seq11, uint32_t& seq_len) {

    uint32_t seq_length = times * len;
    uint8_t* label = (uint8_t*)malloc(sizeof(uint8_t) * seq_length);
    if (!label) {
        cerr << "Memory allocation failed!" << endl;
        exit(EXIT_FAILURE);
    }
    memset(label, 0, seq_length * sizeof(uint8_t));
    random_device rd;
    mt19937 gen(rd());
    int edit_distance = seq_length * edit_distance_radio;
    if (seq_length < edit_distance) {
        free(label); // 释放内存
        return;
    }

    uniform_int_distribution<> dis_label(1, 3);
    uniform_int_distribution<> dis_insert(0, seq_length);
    int32_t edit_length = 0;
    int* index = (int*)malloc(sizeof(int) * seq_length);
    if (!index) {
        cerr << "Memory allocation failed!" << endl;
        free(label);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < seq_length; i++) {
        index[i] = i;
    }

    srand(time(NULL)); // 设置随机种子
    shuffle(index, seq_length);

    for (size_t i = 0; i < edit_distance; i++) {
        label[index[i]] = dis_label(gen); // 随机修改为1, 2, 或3, 1删除，2替换，3插入
        edit_length = edit_length + (label[index[i]] - 2);
    }

    int seq_length1 = start_len + edit_length + seq_length;
//    cout << start_len << ' ' << edit_length << ' ' << seq_length << ' ' << seq_length1 << endl;

    int insert_pos = dis_insert(gen);
    char* seq1 = NULL;
    char* seq = (char*)malloc(sizeof(char) * (seq_length1 + 1));
    if (!seq) {
        cerr << "Memory allocation failed!" << endl;
        free(label);
        free(index);
        exit(EXIT_FAILURE);
    }
    memset(seq, '\0', seq_length1 + 1);
    char* temp = NULL;

    generate_seq(insert_pos, &temp);
    for (int i = 0; i < insert_pos; i++) {
        seq[i] = temp[i];
    }
    free(temp); // 释放temp

    generate_perfect_seq(len, times, &seq1);

    uint32_t i = 0, j = 0;
    for (; i < seq_length; i++) {
        if (label[i] == 0) {
            seq[insert_pos + j] = seq1[i];
        } else if (label[i] == 1) {
            continue; // 删除操作
        } else if (label[i] == 2) {
            generate_seq(1, &temp);
            seq[insert_pos + j] = temp[0];
            free(temp); // 释放temp
        } else if (label[i] == 3) {
            generate_seq(1, &temp);
            seq[insert_pos + j] = temp[0];
            j++;
            seq[insert_pos + j] = seq1[i];
            free(temp); // 释放temp
        }
        j++;
    }

    generate_seq(start_len - insert_pos, &temp);
    for (uint32_t i = 0; i < start_len - insert_pos; i++) {
        seq[insert_pos + j] = temp[i];
        j++;
    }

    *seq11 = seq;
    seq_len = seq_length1;

    // 释放分配的内存
    free(temp);
    free(seq1);
    free(label);
    free(index);
//    free(seq); // 释放seq
}
//筛选出数量较好的radio
vector<int> find_top_radio_index(vector<vector<double>> data, int seq_length){
	int top_n=1;
	vector<pair<double, int>>radios;
	int i=0;
	for(const auto& vec:data){
		if(vec.size()<10) {
			i++;
			continue;
		}
		int count = std::count_if(vec.begin(), vec.end(), [](double value) {
			return value >= 0.95 && value <=1.05; // 判断条件
		});
		double ratio_in_data = static_cast<double>(count) / seq_length; // 计算整数1占据比例
		radios.emplace_back(ratio_in_data, i); // 保存比例和索引
		i++;
	}
	sort(radios.rbegin(),radios.rend());
	vector<int> top_index;
	for(int i=0;i<min(top_n, static_cast<int>(radios.size()));i++){
		top_index.push_back(radios[i].second);
//		if(radios[i].second<<radio_bar) break;
//		radio_1_in= radios[i].first;
//		cout<<radios[i].first*seq_length/data[radios[i].second].size()<<' ';
		ratios_1_in.push_back(radios[i].first*seq_length/data[radios[i].second].size());
	}
	double sum = 0.0;
	int count = 0;
	// 遍历前十个元素
//	for (size_t i = 0; i < radios.size() && i < 10; ++i) {
//		sum += radios[i].first; // 累加 .first 值
//		count++; // 计数
//	}
//	radio_1_in=count > 0 ? sum / count : 0.0;
	return top_index;
}
string generateInfoFilePath(const string fasta_file) {
    string info_file = fasta_file;
    // 查找 ".fasta" 后缀
    size_t pos = info_file.rfind(".fasta");
    if (pos != string::npos && pos == info_file.length() - 6 ) {
          // 如果找到了 ".fasta" 后缀，并且它在末尾
        info_file.erase(pos); // 删除 ".fasta" 后缀
    }
    info_file += ".info"; // 添加 "_info.txt"
    return info_file;
}
string generate_unique_filename(const string& base_path, int len, int time) {
    // 构造初始文件路径
    std::string file_name = base_path + std::to_string(len) + "_" + std::to_string(time) + ".fasta";
    std::ifstream file(file_name);
    int counter = 1;

    // 如果文件已存在，生成新的路径
    while (file.good()) {
        file_name = base_path + std::to_string(len) + "_" + std::to_string(time) + "(" + std::to_string(counter) +")"+ ".fasta";
        file.close();
        file.open(file_name);
        ++counter;
    }

    return file_name;
}
int generateRandom(int lower, int upper) {
    int value= rand() % (upper - lower + 1) + lower;
	// if (value > 101) {
	// 	return (value + 50) / 100 * 100; // 向最近的 100 取整
	// } else 
	if (value > 10) {
		return (value + 5) / 10 * 10; // 向最近的 10 取整
	}
	return value; // 小于或等于 10 的数字保持不变
}
int generateRandom(int lower, int upper, bool flag) {
    int value= rand() % (upper - lower + 1) + lower;
	return value; // 小于或等于 10 的数字保持不变
}

void generate(vector<int> len_lows, vector<int> len_ups, int time_low, int time_up,  int cell_num, int seq_num, 
	double ed_low,double ed_up, char* init_seq_file,vector< bool > have_flags, char* distination_file){
	char* init_seq;
	uint64_t init_seq_len=0;
	int len,time_i;
	double ed;
	vector<int> lens;
	vector<int> times;
	vector<int> insert_pos;
	vector<double> edit_distance;
	int insert_p=0;
	int cell_insert_len=0;
	if(have_flags[0]==0){
		init_seq_len=2*len_ups[0]*time_up;;
		cell_insert_len= init_seq_len/cell_num;
		// insert_p=generateRandom(0, cell_insert_len,1);
		cout<<insert_p<<endl;
	}else {
		ReadSeq_ref(&init_seq, &init_seq_len, init_seq_file);
		cell_insert_len= init_seq_len/cell_num;
		// insert_p=generateRandom(0, cell_insert_len,1);
	}
	if(have_flags[1]==0){
		seq_num=generateRandom(10,20);
	}
	srand(time(0)); // 随机种子

	//生成命名文件
	string base_path="generate_fasta/";
	string file_name_s;
	string info_file_name_s;
	if(distination_file== NULL){
		file_name_s=generate_unique_filename(base_path, len_ups[0], time_up);
		info_file_name_s=generateInfoFilePath(file_name_s);
	}else{
		string temp(distination_file);
		file_name_s= base_path+"fasta/"+temp+".fasta";
		info_file_name_s = base_path+"info/"+temp+".info";
		cout<<file_name_s<<endl;
	}
	const char* file_name=file_name_s.c_str();
	const char* info_file_name = info_file_name_s.c_str();
	cout<<"generate result file name is: \t"<<file_name<<endl;
	cout<<"generate result info_file name is: \t"<<info_file_name<<endl;
	write_file("", file_name, "w");
	write_file("", info_file_name, "w");
	//将生成数据写入文件
	for(int i=1;i<=seq_num;i++){
		insert_p=0;
		lens.clear();
		times.clear();
		edit_distance.clear();
		insert_pos.clear();
		//生成模拟序列数据
		if(have_flags[2]==0){
			cell_num=generateRandom(2,10);//重复单元个数
		}
		for(int j=0;j<cell_num;j++){
			if(have_flags[3]==0){
				len=generateRandom(10,100);
			}else {
				for(int k= 0;k< len_lows.size();k++){
					len=generateRandom(len_lows[k],len_ups[k]);
					lens.push_back(len);
				}
			}
			if(have_flags[4]==0){
				time_i= generateRandom(10,100);
			}else{
				for(int k= 0;k< len_lows.size();k++){
					time_i= generateRandom(time_low, time_up);
					times.push_back(time_i);
				}
			}
			if(have_flags[5]==0){
				ed=generateRandom(0,10)/100.0;
			}else{
				for(int k= 0;k< len_lows.size();k++){
					ed=generateRandom(int(ed_low*100),int(ed_up*100),0)/100.0;
					edit_distance.push_back(ed);
				}
			}
			for(int k= 0;k< len_lows.size();k++){
				insert_p=generateRandom(insert_p, cell_insert_len*(j+1), 1);
				insert_pos.push_back(insert_p);
			}
			// cout<<init_seq_len<<' '<<seq_num<<' '<<cell_num<<' '<<len<<' '<<time_i<<' '<<ed<<endl;
		}
		
		//生成模拟序列
		uint32_t seq_len=0;
		char* seq=NULL;
		generate_seq_list_ns(init_seq_len, insert_pos, lens, times, edit_distance, &seq, seq_len, have_flags[0], init_seq);
		char write_buffer[5000];

		sprintf(write_buffer, "> %d %zu  \n", i, seq_bound.size());
		write_file(write_buffer,file_name,"a");
		// write_file(write_buffer,info_file_name,"a");
		for(int j=0; j<seq_bound.size();j++){
			sprintf(write_buffer, "> %d\t%d\t%d\t%5d\t%5d\t%d\t%s\n", i,j+1, lens[j], seq_bound[j].first, seq_bound[j].second,
				 times[j],  seq_cell[j].c_str());
			write_file(write_buffer,info_file_name,"a");
		}
		write_file(seq, file_name,"a");
		// cout<<seq<<endl;
		sprintf(write_buffer, "\n");
		write_file(write_buffer,file_name,"a");
		seq_cell.clear();
		free(seq);
	}
}
int write_file(char* data, const char* file_name, const char * flag){
	
	fs::path dirpath = fs::path(file_name).parent_path();
	if(!fs::exists(dirpath)){
		fs::create_directories(dirpath);
	}
	FILE *file = fopen(file_name, flag);
	if (file == NULL) {
		perror("无法打开文件");
		return 1;
	}
	if (fputs(data, file) == EOF) {
		perror("写入文件失败");
		fclose(file);
		return 1;
	}
	fclose(file);
	return 1;
}
int write_file(char* data, const char* file_name){
	
	fs::path dirpath = fs::path(file_name).parent_path();
	if(!fs::exists(dirpath)){
		fs::create_directories(dirpath);
	}
	FILE *file = fopen(file_name,"a");
	if (file == NULL) {
		perror("无法打开文件");
		return 1;
	}
	if (fputs(data, file) == EOF) {
		perror("写入文件失败");
		fclose(file);
		return 1;
	}
	fclose(file);
	return 1;
}
//void plot_data_by_py(vector<double> data, char* kmer,vector<uint32_t> index){
//	vector<double> y=data;
//	vector<uint32_t> x(y.size());
//	for(int i=0;i<y.size();i++){
//		x[i]=i;
//	}
//	if (x.empty() || y.empty() || x.size() != y.size()) {
//		std::cerr << "Error: x and y must be non-empty and of the same size." << std::endl;
//		return;
//	}
//
//
//	plt::scatter(x, y, 1, {{"cmap", "viridis"}, {"alpha", "0.5"}});
////	plt::imshow({{0, 1}}, {{"cmap", "viridis"}, {"visible", false}});  // 用作颜色条
//	plt::title("Scatter Plot of Simulated Ratio Features");
//	plt::xlabel("Position");
//	plt::ylabel("Ratio");
//	std::string filename = string(kmer) + "_ratio.png"; // 创建文件名
//	plt::save(filename);
//
//	std::cout << "图形已保存为 " << filename << std::endl;
//}



