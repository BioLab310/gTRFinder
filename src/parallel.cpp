
#include "parallel.h"
void* seq_index_parallel(void* arg){
	thread_data* a= (thread_data*) arg;
//	vector<vector<int>>& index=a->index;
//	ReadSeq_parallel(a->fp,a->seq, a->begin, a->end);
	generate_index_int(a->index, a->seq, a->begin, a->end, a->kmer_len);
	pthread_exit(NULL);
}


void index_gen_parallel(char* pathFile, uint32_t kmer_len,int thread )
{
//	pthread_mutex_init(&mutex_index, NULL);
	if (kmer_len > 15) {
		cout << "error" << endl;
		return ;
	}
	pthread_t* t;
	int thread_num=thread;
	t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
	thread_data* data=(thread_data*)malloc(sizeof(thread_data)*thread_num);
	//1. 读取序列文件
	struct RefFilePath p_ref_path;
	getRefFilePathes(pathFile, p_ref_path);
	cout<<"dataset1 filenum is "<<p_ref_path.NumberOfPathes<<endl;
	for (uint32_t ref_i = 0; ref_i < p_ref_path.NumberOfPathes; ref_i++) {
		//2. 初始化hash table for saving the (k+1)-mers
		char* seq = (char*)malloc(sizeof(char) * 1);
		uint64_t seq_length = 0;
		FILE* fp;
		fp = fopen(p_ref_path.pRefFilePath[ref_i],"r+");
//		uint32_t buffer_size = 256;
//		char buffer_line[256];
//		memset(buffer_line, 0, buffer_size);
//		if (fp == NULL)
//		{
//			cout << "file can not be open!" << endl;
//			return;
//		}
//		while (fgets(buffer_line, buffer_size - 1, fp) != NULL)
//		{
//			if (buffer_line[0] == '>')
//				continue;
//			else
//			{
//				seq_length+=strlen(buffer_line);
//			}
//			memset(buffer_line, 0, buffer_size);
//		}
		struct timeval tvs, tve;
		double time_token=0;
		uint32_t p_size = pow(4, kmer_len);
		vector < vector<uint32_t>>p_int(p_size);
		ReadSeq_ref(&seq, &seq_length, p_ref_path.pRefFilePath[ref_i]);
		uint32_t uint=seq_length/thread_num;
		if(uint==0||uint==1){
			generate_index_int(p_int, seq,0, seq_length, kmer_len);
		}
		else{
			gettimeofday(&tvs, NULL);
			seq=(char*)malloc(sizeof(char)*seq_length+1);
			for(int i=0;i<thread_num;i++){
				uint32_t begin;
				uint32_t end;
				if(i==0) {
					begin=0;
					end=(i+1)*uint;
				}
				else if(i==thread_num-1){
					begin=end-kmer_len+1;
					end=seq_length;
				}
				else{
					begin=end-kmer_len+1;
					end=(i+1)*uint;
				}
//				data[i].ref_path=p_ref_path.pRefFilePath[ref_i];
//				data[i].index= p_int;
//				data[i].seq=seq;
//				data[i].begin= begin;
//				data[i].end=end;
//				data[i].kmer_len=kmer_len;
				new(data+i) thread_data(fp,p_int, seq, begin, end, kmer_len);
				if(pthread_create(t+i, NULL, seq_index_parallel, (void*)(data+i))!=0){
					cout<< "error"<<endl;
				}
			}
			for(int i=0;i< thread_num;i++){
				pthread_join(t[i],NULL);
			}
			//合并数据
			for(int p_i=0;p_i<thread_num;p_i++){
				for(int i=0;i<data[p_i].index.size();i++){
					for(int j=0;j<data[p_i].index[i].size();j++){
						p_int[i].push_back(data[p_i].index[i][j]);
					}
				}
			}
		}
		gettimeofday(&tve, NULL);
//		time_token= (tve.tv_sec-tvs.tv_sec)*1e6;
//		time_token= (time_token+ tve.tv_usec- tvs.tv_usec)*1e-6;
//		cout<< "parallel_index_trans take time: "<<time_token<<endl;
//		index_display_int(p_int, kmer_len);
//		cout<<seq<<endl;
		free(seq);
		free(data);
		free(t);
		vector<vector<double>> charactor;
		character_cmp(p_int, charactor);
	}
	freeRefFilePath(&p_ref_path);
}
//降低函数间的耦合度
void index_gen_parallel_seq1(char* seq, uint32_t seq_length, uint32_t kmer_len,int thread, vector<vector<uint32_t>>& p_int)
{
//	pthread_mutex_init(&mutex_index, NULL);
	if (kmer_len > 15) {
		cout << "error" << endl;
		return ;
	}
	pthread_t* t;
	int thread_num=thread;
	t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
	thread_data* data = new thread_data[thread_num];
	struct timeval tvs, tve;
	double time_token=0;
	uint32_t p_size = pow(4, kmer_len);
	p_int.resize(p_size);
	uint32_t uint=seq_length/thread_num;
	if(uint==0||uint==1){
		generate_index_int(p_int, seq,0, seq_length, kmer_len);
	}
	else{
		
		gettimeofday(&tvs, NULL);
		for(int i=0;i<thread_num;i++){
			uint32_t begin;
			uint32_t end;
			if(i==0) {
				begin=0;
				end=(i+1)*uint;
			}
			else if(i==thread_num-1){
				begin=end-kmer_len+1;
				end=seq_length;
			}
			else{
				begin=end-kmer_len+1;
				end=(i+1)*uint;
			}
			new(data+i) thread_data(NULL,p_int, seq, begin, end, kmer_len);
			if(pthread_create(t+i, NULL, seq_index_parallel, (void*)(data+i))!=0){
				cout<< "error"<<endl;
			}
		}
		for(int i=0;i< thread_num;i++){
			pthread_join(t[i],NULL);
		}
		//合并数据
		for(int p_i=0;p_i<thread_num;p_i++){
			for(int i=0;i<data[p_i].index.size();i++){
				for(int j=0;j<data[p_i].index[i].size();j++){
					p_int[i].push_back(data[p_i].index[i][j]);
				}
			}
		}
	}
	gettimeofday(&tve, NULL);
//	vector<int> top_radios=find_top_radio_index(charactor,seq_length);//返回的是列表序号
//	true_false_compute(p_int[top_radios[0]]);


//	index_display_int(p_int,kmer_len,charactor, seq_length);
	delete[] data;
	free(t);
}
void index_gen_parallel_seq(char* seq, uint32_t seq_length, uint32_t kmer_len,int thread)
{
//	pthread_mutex_init(&mutex_index, NULL);
	if (kmer_len > 15) {
		cout << "error" << endl;
		return ;
	}
	pthread_t* t;
	int thread_num=thread;
	t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
	thread_data* data = new thread_data[thread_num];
	struct timeval tvs, tve;
	double time_token=0;
	uint32_t p_size = pow(4, kmer_len);
	vector < vector<uint32_t>>p_int(p_size);
	uint32_t uint=seq_length/thread_num;
	if(uint==0||uint==1){
		generate_index_int(p_int, seq,0, seq_length, kmer_len);
	}
	else{
		gettimeofday(&tvs, NULL);
		for(int i=0;i<thread_num;i++){
			uint32_t begin;
			uint32_t end;
			if(i==0) {
				begin=0;
				end=(i+1)*uint;
			}
			else if(i==thread_num-1){
				begin=end-kmer_len+1;
				end=seq_length;
			}
			else{
				begin=end-kmer_len+1;
				end=(i+1)*uint;
			}
			new(data+i) thread_data(NULL,p_int, seq, begin, end, kmer_len);
			if(pthread_create(t+i, NULL, seq_index_parallel, (void*)(data+i))!=0){
				cout<< "error"<<endl;
			}
		}
		for(int i=0;i< thread_num;i++){
			pthread_join(t[i],NULL);
		}
		//合并数据
		for(int p_i=0;p_i<thread_num;p_i++){
			for(int i=0;i<data[p_i].index.size();i++){
				for(int j=0;j<data[p_i].index[i].size();j++){
					p_int[i].push_back(data[p_i].index[i][j]);
				}
			}
		}
	}
	gettimeofday(&tve, NULL);
	vector<vector<double>> charactor;
	character_cmp(p_int, charactor);
	vector<int> top_radios=find_top_radio_index(charactor,seq_length);//返回的是列表序号
	true_false_compute(p_int[top_radios[0]]);


//	index_display_int(p_int,kmer_len,charactor, seq_length);
	delete[] data;
	free(t);
}
void index_gen_parallel_seq_plot(char* seq, uint32_t seq_length, uint32_t kmer_len,int thread, int flag)
{
//	pthread_mutex_init(&mutex_index, NULL);
	if (kmer_len > 15) {
		cout << "error" << endl;
		return ;
	}
	pthread_t* t;
	int thread_num=thread;
	t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);
	thread_data* data = new thread_data[thread_num];
	struct timeval tvs, tve;
	double time_token=0;
	uint32_t p_size = pow(4, kmer_len);
	vector < vector<uint32_t>>p_int(p_size);
	uint32_t uint=seq_length/thread_num;
	if(uint==0||uint==1){
		generate_index_int(p_int, seq,0, seq_length, kmer_len);
	}
	else{
		gettimeofday(&tvs, NULL);
		for(int i=0;i<thread_num;i++){
			uint32_t begin;
			uint32_t end;
			if(i==0) {
				begin=0;
				end=(i+1)*uint;
			}
			else if(i==thread_num-1){
				begin=end-kmer_len+1;
				end=seq_length;
			}
			else{
				begin=end-kmer_len+1;
				end=(i+1)*uint;
			}
			new(data+i) thread_data(NULL,p_int, seq, begin, end, kmer_len);
			if(pthread_create(t+i, NULL, seq_index_parallel, (void*)(data+i))!=0){
				cout<< "error"<<endl;
			}
		}
		for(int i=0;i< thread_num;i++){
			pthread_join(t[i],NULL);
		}
		//合并数据
		for(int p_i=0;p_i<thread_num;p_i++){
			for(int i=0;i<data[p_i].index.size();i++){
				for(int j=0;j<data[p_i].index[i].size();j++){
					p_int[i].push_back(data[p_i].index[i][j]);
				}
			}
		}
	}
	gettimeofday(&tve, NULL);
	vector<vector<double>> charactor;
	character_cmp(p_int, charactor);
	index_display_int(p_int,kmer_len,charactor, seq_length, flag);
	delete[] data;
	free(t);
}
void ReadSeq_parallel(FILE* fp,char* seq, uint64_t begin, uint64_t end){
//	FILE* fp;
//	cout<< p_ref_path<<endl;
//	fp= fopen(p_ref_path, "r+");
	if(fp==NULL){
		cout<<"file is not openned"<<endl;
		return;
	}
	fseek(fp,0,SEEK_SET);
	char ch;
	while((ch=fgetc(fp)!='\n'));
	fseek(fp,0,SEEK_CUR);
	fseek(fp,begin, SEEK_CUR);
	uint64_t len=end-begin;
	fread(seq+begin, sizeof(char), len,fp);
}

