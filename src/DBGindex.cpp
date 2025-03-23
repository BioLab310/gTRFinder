#include "DBGindex.h"
uint64_t M_nrp_lenl = 1024 * 1024;
uint32_t add_len = 1024;
void index_gen_disp(char* pathFile, uint32_t kmer_len)
{
	if (kmer_len > 11) {
		cout << "error" << endl;
	}
	//1. 读取序列文件
	struct RefFilePath p_ref_path;
	getRefFilePathes(pathFile, p_ref_path);
	for (uint32_t ref_i = 0; ref_i < p_ref_path.NumberOfPathes; ref_i++) {
		//2. 初始化hash table for saving the (k+1)-mers
		char* seq = (char*)malloc(sizeof(char) * 1);
		uint64_t seq_length = 0;
//		struct bit64Hash* p_root = NULL;
//		p_root = initial256BitHashTable();
//		generate_index(p_ref_path.pRefFilePath[ref_i], p_root, seq, seq_length, kmer_len);
//		index_display(p_root, kmer_len);
//		delete64BitHash(p_root, HashSize);
		if (kmer_len > 15) return;
		uint32_t p_size = pow(4, kmer_len);
		vector < vector<uint32_t>>p_int(p_size);
		ReadSeq_ref(&seq, &seq_length, p_ref_path.pRefFilePath[ref_i]);
		struct timeval tvs, tve;
		double time_token=0;
		gettimeofday(&tvs, NULL);
		generate_index_int(p_int, seq,0, seq_length, kmer_len);
		gettimeofday(&tve, NULL);
		time_token= (tve.tv_sec-tvs.tv_sec)*1e6;
		time_token= (time_token+ tve.tv_usec- tvs.tv_usec)*1e-6;
		cout<< "no_parallel_index_trans take time: "<<time_token<<endl;
//		index_display_int(p_int, kmer_len);
//		cout<<seq<<endl;
		vector<vector<double>> character;
		character_cmp(p_int, character);
//		index_display_int(p_int,kmer_len, character);
		free(seq);

	}
	freeRefFilePath(&p_ref_path);
}
void index_gen_seq(char* seq, uint32_t seq_length,uint32_t kmer_len)
{
	if (kmer_len > 11) {
		cout << "error" << endl;
	}
	if (kmer_len > 15) return;
	uint32_t p_size = pow(4, kmer_len);
	vector < vector<uint32_t>>p_int(p_size);
	struct timeval tvs, tve;
	double time_token=0;
	gettimeofday(&tvs, NULL);
	generate_index_int(p_int, seq,0, seq_length, kmer_len);
	gettimeofday(&tve, NULL);
	time_token= (tve.tv_sec-tvs.tv_sec)*1e6;
	time_token= (time_token+ tve.tv_usec- tvs.tv_usec)*1e-6;
	cout<< "no_parallel_index_trans take time: "<<time_token<<endl;
	vector<vector<double>> character;
	character_cmp(p_int, character);
//	index_display_int(p_int,kmer_len, character);
//	free(seq);
}
void generate_index(char* p_ref_path, bit64Hash* p_root, char* seq, uint64_t seq_length, uint32_t kmer_len) {
	uint64_t kmer_int = 0;
	cal_hash_value_directly_64bit(seq, kmer_int, kmer_len);
	insert64BitHashTable(p_root, kmer_int, 0);
	for (int i = 1; i < seq_length - kmer_len; i++) {
		cal_hash_value_indirectly_64bit(seq+i, kmer_int, kmer_int, kmer_len);
		insert64BitHashTable(p_root, kmer_int, i);
	}
}
void index_display(bit64Hash* ph,uint32_t kmer_len) {
	char* kmer = (char*)malloc(sizeof(char) * kmer_len);
	for (int i = 0; i < HashSize; i++) {
		if (ph->len[i] == 0) continue;
		for (int j = 0; j < ph->len[i]; j++) {
			int_value_to_kmer(ph->p_kmer[i][j], kmer, kmer_len);
			cout << kmer<<" of index is: ";
			for (int k = 0; k < ph->p_pos[i][j].size(); k++) {
				cout << ph->p_pos[i][j][k] << ' ';
			}//索引值
			cout << endl;
		}
	}
	free(kmer);
}

void generate_index_int(vector<vector<uint32_t>>& p,char* seq, uint64_t begin,uint64_t end, uint32_t kmer_len) {
//	ReadSeq_ref(&seq, &seq_length, p_ref_path);
	uint64_t kmer_int = 0;
	cal_hash_value_directly_64bit(seq+begin, kmer_int, kmer_len);
	p[kmer_int].push_back(begin);
	for (int i = begin+1; i < end- kmer_len+1; i++) {
		cal_hash_value_indirectly_64bit(seq+i, kmer_int, kmer_int, kmer_len);//窗口的前一个指针向后移动
		p[kmer_int].push_back(i);
	}
}
void index_display_int(vector<vector<uint32_t>> p, uint32_t kmer_len,vector<vector<double>> ratio, int seq_length, int flag) {
	char* kmer = (char*)malloc(sizeof(char) * kmer_len);
	vector<int> top_radios=find_top_radio_index(ratio,seq_length);//返回的是列表序号
	int ii=0;
	for(auto i:top_radios){
		if (p[i].size() == 0) continue;
		int_value_to_kmer(i, kmer, kmer_len);
		plot_data_by_py(ratio[i],kmer_radio,ii, p[i], flag);
		ii++;
	}
	free(kmer);
}

void character_cmp(vector<vector<uint32_t>> index,vector<vector<double>>& character){
	//初始化数据
	vector<vector<int>> diff;
	diff.resize(index.size());
	character.resize(index.size());
	for(int i=0;i<index.size();i++){
		if(index[i].size()<2) continue;
		diff[i].resize(index[i].size()-1);
		if(index[i].size()<3) continue;
		character[i].resize(index[i].size()-2);
//		cout<<"index: "<<index[i].size()<<endl;
	}
	for(int i=0;i<index.size();i++){
		for(int j=0; j<diff[i].size();j++){
			diff[i][j]=index[i][j+1]-index[i][j];
		}
	}
	for(int i=0;i<index.size();i++){
		for(int j=0; j<character[i].size();j++){
			character[i][j]=diff[i][j+1]*1.0/diff[i][j];
		}
	}
//	return compute_space(diff, character);
}

//void generate_index(char* p_ref_path ,bit64Hash* p_root, char* nrp, uint64_t& nrp_lenl, char* seq, uint64_t seq_length,uint32_t kmer_len) {
//	uint64_t seq_current = 0;
//	ReadSeq_ref(&seq, &seq_length, p_ref_path);
//	//读取每个文件中的acgt序列
//	char* tmp_seq;
//	tmp_seq = seq;
//	cal_hash_value_directly_64bit(tmp_seq, seq_current, kmer_len);
//	//生成第一个k_mer
//	int32_t arrayID = search64BitHashTable(p_root, seq_current);
//	if (arrayID == -1 || nrp_lenl == 0) {
//		insert64BitHashTable(p_root, seq_current, 0);
//		if (nrp[nrp_lenl - 1] == '$' || nrp_lenl == 0)
//		{
//			if (nrp_lenl + kmer_len + 1 >= M_nrp_lenl)
//			{
//				nrp = (char*)realloc(nrp, sizeof(char) * (M_nrp_lenl + add_len));
//				M_nrp_lenl = M_nrp_lenl + add_len;
//
//			}
//			for (uint32_t i = 0; i < kmer_len + 1; i++)
//			{
//				nrp[nrp_lenl] = tmp_seq[i];
//				nrp_lenl++;
//			}
//		}
//		else
//		{
//			// in this case nrp[nrp_lenl-1]!='$'&& nrp_lenl!=0 && the first (k+1)-mer of a new reference
//			if (nrp_lenl + 1 + kmer_len + 1 >= M_nrp_lenl)
//			{
//				nrp = (char*)realloc(nrp, sizeof(char) * (M_nrp_lenl + add_len));
//				M_nrp_lenl = M_nrp_lenl + add_len;
//
//			}
//			nrp[nrp_lenl] = '$';
//			nrp_lenl++;
//			for (uint32_t i = 0; i < kmer_len + 1; i++)
//			{
//				nrp[nrp_lenl] = tmp_seq[i];
//				nrp_lenl++;
//			}
//			// the following two line are the original operation
////				nrp[nrp_lenl]=temp_seq[kmer_len];
////				nrp_lenl++;
//		}
//	}
//	else
//	{
//		insert_pos_exist(p_root, arrayID, seq_current, 0);
//		if (nrp[nrp_lenl - 1] != '$')
//		{
//			nrp[nrp_lenl] = '$';
//			nrp_lenl++;
//		}
//	}
//
//	//后续k_mer
//	for (uint32_t ii = 1; ii < seq_length - kmer_len+1; ii++) {
//		cal_hash_value_indirectly_64bit(tmp_seq + ii, seq_current, seq_current, kmer_len);
//		int32_t arrayID = search64BitHashTable(p_root, seq_current);
//		if (arrayID == -1) {
//			insert64BitHashTable(p_root, seq_current, ii );
//			if (nrp[nrp_lenl - 1] == '$' || nrp_lenl == 0)
//			{
//				if (nrp_lenl + kmer_len + 1 >= M_nrp_lenl)
//				{
//					nrp = (char*)realloc(nrp, sizeof(char) * (M_nrp_lenl + add_len));
//					M_nrp_lenl = M_nrp_lenl + add_len;
//
//				}
//				for (uint32_t i = 0; i < kmer_len + 1; i++)
//				{
//					nrp[nrp_lenl] = tmp_seq[ii + i];
//					nrp_lenl++;
//				}
//			}
//			else
//			{
//				if (nrp_lenl + 1 >= M_nrp_lenl)
//				{
//					nrp = (char*)realloc(nrp, sizeof(char) * (M_nrp_lenl + add_len));
//					M_nrp_lenl = M_nrp_lenl + add_len;
//
//				}
//				nrp[nrp_lenl] = tmp_seq[ii + kmer_len];
//				nrp_lenl++;
//			}
//		}
//		else {
//			insert_pos_exist(p_root, arrayID, seq_current, ii );
//			if (nrp[nrp_lenl - 1] != '$')
//			{
//				nrp[nrp_lenl] = '$';
//				nrp_lenl++;
//			}
//		}
//	}
//
//
//}
//void index_display(DBGindex *dbg,uint32_t kmer_len) {
//	for (int i = 0; i < dbg->len; i++) {
//		for (int j = 0; j < kmer_len; j++) {
//			cout << dbg->kmer[i][j] ;
//		}
//		cout << " index is ";
//		for (int j = 0; j < dbg->index[i].size(); j++) {
//			cout << dbg->index[i][j] << ' ';
//		}
//		cout << endl;
//	}
//}
//void store_index(char* nrp, bit64Hash* ph ,uint32_t k_mer_len,uint32_t nrp_lenl, DBGindex* dbg) {
//	const char* delimiter = "$";
//	char* context = nullptr;
//	char* kmer = (char*)malloc(sizeof(char) * k_mer_len);
//	//最后一个数据为空
//	uint32_t len= 0;
//	for (int i = 0; i < nrp_lenl - k_mer_len; i++) {
//		if (nrp[i + 3] == '$') {
//			if (i + 4 >= nrp_lenl - k_mer_len) break;
//			i = i + 4;
//		}
//		for (int j = i; j < i + k_mer_len; j++) {
//			kmer[j - i] = nrp[j];
//		}
//		len++;
//	}
//	dbg->kmer = (char**)malloc(sizeof(char*) * len);
//	for (int i = 0; i < len; i++) {
//		dbg->kmer[i] = (char*)malloc(sizeof(char) * k_mer_len);
//	}
//	dbg->index.resize(len);
//	dbg->len = len;
//	len = 0;
//	for(int i = 0; i < nrp_lenl-k_mer_len;i++){
//		if (nrp[i + 3] == '$') {
//			if (i + 4 >= nrp_lenl - k_mer_len) break;
//			i = i + 4;
//		}
//		for (int j = i; j < i+k_mer_len; j++) {
//			kmer[j - i] = nrp[j];
//		}
//		uint64_t current = 0;
//		cal_hash_value_directly_64bit(kmer, current, k_mer_len);
//		uint64_t shift = bit64HashFunction(current);//找到hash位置
//		int32_t arrayID = search64BitHashTable(ph, current);
//		for (int j = 0; j < k_mer_len; j++) {
//			dbg->kmer[len][j] = kmer[j];
//		}
//		dbg->index[len] = ph->p_pos[shift][arrayID];
//		len++;
//	}
//}
//void index_display(char* nrp, bit64Hash* ph, uint32_t k_mer_len, uint32_t nrp_lenl) {
//	const char* delimiter = "$";
//	char* context = nullptr;
//	char* kmer = (char*)malloc(sizeof(char) * k_mer_len);
//	//最后一个数据为空
//	for (int i = 0; i < nrp_lenl - k_mer_len; i++) {
//		if (nrp[i + 3] == '$') {
//			if (i + 4 >= nrp_lenl - k_mer_len) break;
//			i = i + 4;
//		}//出现重复段，跳过
//		for (int j = i; j < i + k_mer_len; j++) {
//			kmer[j - i] = nrp[j];
//		}
//		uint64_t current = 0;
//		cal_hash_value_directly_64bit(kmer, current, k_mer_len);
//		uint64_t shift = bit64HashFunction(current);//找到hash位置
//		int32_t arrayID = search64BitHashTable(ph, current);
//		cout << kmer << ' ' << "index is: ";
//		for (int ii = 0; ii < ph->p_pos[shift][arrayID].size(); ii++) {
//			cout << ph->p_pos[shift][arrayID][ii] << ' ';
//		}
//		cout << endl;
//
//	}
//}
