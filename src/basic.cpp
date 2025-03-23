

#include "basic.h"
pthread_mutex_t mutex_index;
vector<pair<int,int>> seq_bound;
int len_global=0;
int time_global=0;
double edit_distance_radio_global=0;
int init_seq_time=0;
vector<double> ratios_1_in;
vector<double> false_charactor_ratio;
string file_name_s;
Buffer buffer;
void initBuffer(Buffer& buffer) {
    buffer.seq = nullptr;
    buffer.length = 0;
    buffer.pos.clear();
}


void get_para(struct bit256KmerPara* para1, uint32_t kmer_length)
{
	para1->kmer1Len = kmer_length * 2;//
	para1->remainer1to64 = para1->kmer1Len % 64;
	para1->kmer64Len = para1->kmer1Len / 64 + (para1->remainer1to64 ? 1 : 0);
	if (para1->remainer1to64 == 0)
	{
		para1->codefor1 = 0;
	}
	else
	{
		para1->codefor1 = 0;
		for (uint32_t i = 0; i < para1->remainer1to64 - 1; i++)
		{
			para1->codefor1 = para1->codefor1 | 1;
			para1->codefor1 = para1->codefor1 << 1;
		}
		para1->codefor1 = para1->codefor1 | 1;
	}
}

void ReadSeq(char** seq1, uint32_t* seq_length, const char* p_ref) {
	uint32_t cursize;
	uint32_t maxsize;
	uint32_t addsize;
	maxsize = pow(2, 15);//数据量可以小点
	addsize = pow(2, 15);
	cursize = maxsize;
	char* seq;
	seq = (char*)malloc(sizeof(char) * maxsize);

	uint32_t buffer_size = 256;
	char buffer_line[256] = { 0 };//the last num is '\0'
	FILE* fp ;
	fp = fopen(p_ref,"r+");
	if (fp == NULL) {
		cout << "file can not be open" << endl;
		return;
	}

	uint32_t len = 0;
	while (fgets(buffer_line, buffer_size, fp) != NULL) {
		if (buffer_line[0] == '>') {
			continue;
		}
		else {
			if (len + buffer_size < cursize) {
				cursize = cursize+0;
			}//cursize is enough ,no change
			else {
				seq = (char*)realloc(seq, sizeof(char) * (cursize + addsize));
				cursize += addsize;
			}//cursize less than the num, expand the space


			for (uint32_t i = 0; i < buffer_size; i++) {
				if (buffer_line[i] == '\n' || buffer_line[i] == '\0') {
					break;
				}//fgets can't read '\n'
				if (buffer_line[i] >= 'a') {
					buffer_line[i] -= 32;
				}
				if (buffer_line[i] != 'A' && buffer_line[i] != 'C' && buffer_line[i] != 'G' && buffer_line[i] != 'T') {
					continue;
				}
				seq[len] = buffer_line[i];
				len++;
			}
		}
		memset(buffer_line,0, buffer_size * sizeof(char));
		*seq_length = len;
		*seq1 = seq;
	}
}

//打开文件路径
void ReadSeq_ref(char** seq1, uint64_t* seq_length, char* p_ref)
{
//	cout<< p_ref<<endl;
	uint32_t buffer_size = 256;
	char buffer_line[256];
	memset(buffer_line, 0, buffer_size);

	FILE* fp;
	fp = fopen(p_ref,"r+");
	if (fp == NULL)
	{
		cout << "file can not be open!" << endl;
		return;
	}
	char ch;
	uint64_t total_size = 0;
	fseek(fp, 0, 2);
	total_size = ftell(fp);//获得文件的字节数
	char* seq;
	seq = (char*)malloc(sizeof(char) * total_size);

	fseek(fp, 0, 0);//指针回到开始

	uint64_t len = 0;
	while (fgets(buffer_line, buffer_size - 1, fp) != NULL)
	{
		if (buffer_line[0] == '>')
			continue;
		else
		{
			for (uint32_t i = 0; i < strlen(buffer_line)&& len<total_size; i++)
			{
				if (buffer_line[i] == '\n' || buffer_line[i] == '\0')
				{
					break;
				}
				if (buffer_line[i] >= 'a')
				{
					buffer_line[i] -= 32;
				}
				if (buffer_line[i] != 'A' && buffer_line[i] != 'C' && buffer_line[i] != 'G' && buffer_line[i] != 'T')
				{
					buffer_line[i] = 'A';
				}
				seq[len] = buffer_line[i];
				len++;
			}
		}
		memset(buffer_line, 0, buffer_size);
	}
	*seq_length = len;
	*seq1 = seq;
//	cout << seq<<endl;
//	cout << "the length of original seq is: " << len << endl;
}
void ReadSeq_ref(string& seq, string& filename) {
    std::ifstream fp(filename.c_str());  // 使用 c_str() 将 string 转换为 const char* 用于 fopen
    if (!fp.is_open()) {
        std::cerr << "file can not be open!" << std::endl;
        return;
    }

    std::stringstream ss;
    std::string line;

    while (std::getline(fp, line)) {
        if (line[0] == '>') {
            continue; // Skip header lines
        }
        else {
            for (char& c : line) {
                if (c == '\n' || c == '\0') {
                    break;
                }
                if (c >= 'a' && c <= 'z') {
                    c = toupper(c); // Convert to uppercase
                }
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
                    c = 'A';  // Replace invalid characters
                }
                ss << c; // Append to the stringstream
            }
        }
    }

    seq = ss.str(); // Assign the content of the stringstream to the output string
    fp.close();
}

//每次计算一个k_mer 的编号,k_mer的长度不大于32
void cal_hash_value_directly_64bit(char* seq, uint64_t &current, \
	uint32_t k_mer_len){
	uint64_t tmp;
	char* k_mer_temp = seq;
	char* loop_tmp = k_mer_temp;
	tmp = 0;
	for (uint32_t j = 0; j < k_mer_len- 1; j++) {
		switch (loop_tmp[j]) {
		case 'A':
			tmp = tmp << 2;
			break;
		case 'C':
			tmp = tmp | 1;
			tmp = tmp << 2;
			break;
		case 'G':
			tmp = tmp | 2;
			tmp = tmp << 2;
			break;
		case 'T':
			tmp = tmp | 3;
			tmp = tmp << 2;
			break;
		default:
			tmp = tmp << 2;
			break;
		}
	}
	switch (loop_tmp[k_mer_len - 1])
	{
		case 'A':
			break;
		case 'C':
			tmp = tmp | 1;
			break;
		case 'G':
			tmp = tmp | 2;
			break;
		case 'T':
			tmp = tmp | 3;
			break;
		default:
			break;
	}
	current = tmp;
}
uint32_t cmp64BitKmer(uint64_t p, uint64_t b) {
	if (p < b) return 0;
	if (p > b) return 1;
	return 2;
}
void cal_hash_value_indirectly_64bit(char* seq,uint64_t &current,uint64_t original,uint32_t k_mer_len){
	char* k_mer_temp = seq;
	uint64_t tmp = 0;
	uint64_t codefor1 = (1ULL << (k_mer_len * 2)) - 1;
	current = original << 2;
	current = current & codefor1;
	switch (k_mer_temp[k_mer_len - 1]) {
	case 'A':
		break;
	case 'C':
		current = current | 1;
		break;
	case 'G':
		current = current | 2;
		break;
	case 'T':
		current = current | 3;
		break;
	default:
		break;
	}
}

uint32_t cmp256BitKmer(uint64_t* a, uint64_t* b, uint32_t len)
{
	uint32_t r=2;
	for(uint32_t i=0;i<len;i++)
	{
		if(a[i]<b[i])
		{
			r=0;
			break;
		}
		else if(a[i]>b[i])
		{
			r=1;
			break;
		}
		else
		{
			if(i==len-1)
			{
				r=2;
				break;
			}
		}
	}
	return r;
}

void getRefFilePathes(char* pathFile, struct RefFilePath& p) {
	ifstream input;
	input.open(pathFile, ios::in);
	uint32_t total_file_num = 0;
	char hv_tmp[128];
	while (input.getline(hv_tmp, 127)) {
		total_file_num++;
	}
	input.close();

	p.NumberOfPathes = total_file_num;
	p.pRefFilePath = (char**)malloc(sizeof(char*) * total_file_num);
	input.open(pathFile, ios::in);
	for (uint32_t i = 0; i < total_file_num; i++) {
		p.pRefFilePath[i] = (char*)malloc(sizeof(char) * 128);
		input.getline(hv_tmp , 127);
		strcpy(p.pRefFilePath[i], hv_tmp);
	}
	input.close();

}
void int_value_to_kmer(uint64_t kmer_int, char* kmer, uint32_t kmer_len) {
	uint64_t mask = 3;
	for (int i = 0; i < kmer_len; i++) {
		uint64_t mask_value = mask << (i * 2);
		uint64_t current = mask_value & kmer_int;
		current = current >> (i * 2);
		switch (current) {
		case 0:
			kmer[kmer_len-i-1] = 'A';
			break;
		case 1:
			kmer[kmer_len - i - 1] = 'C';
			break;
		case 2:
			kmer[kmer_len - i - 1] = 'G';
			break;
		case 3:
			kmer[kmer_len - i - 1] = 'T';
			break;
		}
	}
}
void freeRefFilePath(struct RefFilePath* ref) {
    if (ref->pRefFilePath != NULL) {
        for (uint32_t i = 0; i < ref->NumberOfPathes; ++i) {
            if (ref->pRefFilePath[i] != NULL) {
                free(ref->pRefFilePath[i]);  // 释放每个字符串
            }
        }
        // 释放指针数组本身
        free(ref->pRefFilePath);
    }
//    free(ref);
}

//void init_kti(kmer_TRs_interval& k, int edge_len){
//	k.kmer=(char*)malloc(sizeof(edge_len));
//	k.start.resize(0);
//	k.end.resize(0);
//}

