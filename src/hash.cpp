

#include "hash.h"

struct bit64Hash* initial256BitHashTable()
{
	struct bit64Hash* ph =new bit64Hash;
	ph->p_kmer = (uint64_t**)malloc(sizeof(uint64_t*) * HashSize);
	vector<vector<vector<uint32_t>>> p_pos;
	p_pos.resize(HashSize);
	ph->p_pos = p_pos;

	// 为 len 和 maxlen 分配内存，并将其初始化为 0
	ph->len = (uint64_t*)malloc(sizeof(uint64_t) * HashSize);
	memset(ph->len, 0, sizeof(uint64_t) * HashSize);
	ph->maxlen = (uint64_t*)malloc(sizeof(uint64_t) * HashSize);
	memset(ph->maxlen, 0, sizeof(uint64_t) * HashSize);

	return ph;
}

uint32_t bit64HashFunction(uint64_t current)
{
	uint64_t tmp = 0;
	tmp = current;
	return tmp & (HashSize - 1);//最后20位作为哈希值，k_mer=10
}

int32_t search64BitHashTable(bit64Hash* ph, uint64_t seq_current)
{
	uint64_t shift = bit64HashFunction(seq_current);
	bool label = 1;
	if (ph->len[shift] < MinSizeForSearch) {
		for(uint32_t i = 0; i < ph->len[shift]; i++) {
			if (cmp64BitKmer(ph->p_kmer[shift][i], seq_current) == 2) {
				return i;//返回相同时位置
			}
		}
	}
	else {
		return binary_64_bit_search(ph->p_pos[shift], ph->p_kmer[shift], ph->len[shift], seq_current);
	}
	if(label) return -1;
	return 0;
}


int32_t binary_64_bit_search(vector<vector<uint32_t>> pos, uint64_t* p, uint64_t len, uint64_t key) {
	if (len == 0)
		return -1;
	if (cmp64BitKmer(p[0], key) == 2)
		return 0;
	else if (cmp64BitKmer(p[len - 1], key) == 2)
		return len-1;
	else {
		uint32_t left, mid, right;
		left = 0;
		right = len - 1;
		while (left > right) {
			mid = (left + right) / 2;
			short x = cmp64BitKmer(p[mid], key);
			if (x== 2)
				return mid;
			else if (x == 0)
				left = mid + 1;
			else right = mid - 1;
		}
		return -1;
	}
}
uint32_t pos_insert_64_bit_table(uint64_t* p, uint64_t len, uint64_t k_mer, bool& new_flag) {
	uint32_t r = 0;
	new_flag = true;
	if (len == 0) {
		//new_flag = false;
		return 0;
	}
	if (cmp64BitKmer(p[0], k_mer) == 1 )
		return 0;
	else if (cmp64BitKmer(p[0], k_mer) == 2) {
		new_flag = false;
		return 0;
	}
	else if (cmp64BitKmer(p[len - 1], k_mer) == 0  )
		return len ;
	else if (cmp64BitKmer(p[len - 1], k_mer) == 2) {
		new_flag = false;
		return len;
	}
	else {
		if (len < MinSizeForSearch) {
			for (uint32_t i = 0; i < len; i++) {
				short x = cmp64BitKmer(k_mer, p[i]);
				if (x == 0) {
					return i;
				}
				if (x == 2 ) {
					new_flag = false;
					return i;
				}
			}
		}
		else {
			uint32_t left, mid, right;
			left = 0;
			right = len - 1;
			while (left > right) {
				mid = (left + right) / 2;
				short x = cmp64BitKmer(p[mid],k_mer);
				if (x == 2) {
					new_flag = false;
					return mid;
				}
				else if (x == 0) {
					left = mid + 1;
				}
				else
					right = mid - 1;
			}
			r=left;
		}
	}
	return r;
}
void mvtoNext(uint64_t* p_kmer, vector<vector<uint32_t>>& p_pos,uint32_t pos,uint32_t len) {
	//uint32_t fix = len - 1 + pos;
	for (uint32_t i = pos; i < len; i++) {
		p_kmer[i + 1] = p_kmer[i];
		p_pos[i + 1] = p_pos[i];
	}
}

void insert64BitHashTable(bit64Hash* ph, uint64_t kmer_int, uint32_t pos) {
	uint32_t shift = bit64HashFunction(kmer_int);
	if (ph->maxlen[shift] == 0) {
		ph->p_kmer[shift] = (uint64_t*)malloc(sizeof(uint64_t) * \
			(ph->maxlen[shift] + AddSize));
		ph->p_pos[shift].resize(1);
		ph->maxlen[shift] += AddSize;
	}
	if (ph->len[shift] >= ph->maxlen[shift]) {
		ph->p_kmer[shift] = (uint64_t*)realloc(ph->p_kmer[shift], sizeof(uint64_t) * \
			(ph->maxlen[shift] + AddSize));
		ph->p_pos[shift].resize(ph->maxlen[shift]+AddSize);
		ph->maxlen[shift] += AddSize;
	}
	bool new_flag = true;
	uint32_t insert_pos = pos_insert_64_bit_table(ph->p_kmer[shift], ph->len[shift], kmer_int, new_flag);
	if (new_flag) {
		if (insert_pos != ph->len[shift]) {
			mvtoNext(ph->p_kmer[shift], ph->p_pos[shift], insert_pos, ph->len[shift]);
		}
		ph->p_kmer[shift][insert_pos] = kmer_int;
		vector<uint32_t> p;
		p.push_back(pos);
		ph->p_pos[shift][insert_pos] = p;
		ph->len[shift]++;
	}
	else {
		ph->p_pos[shift][insert_pos].push_back(pos);
	}
}

void deleteBit64Hash(struct bit64Hash* hash, size_t kmer_size) {
    if (hash->p_kmer != nullptr) {
        for (size_t i = 0; i < kmer_size; ++i) {
            delete[] hash->p_kmer[i];  // 释放每个k-mer数组
        }
        delete[] hash->p_kmer;  // 释放p_kmer本身
    }
    delete[] hash->len;
    delete[] hash->maxlen;
    // 释放结构体本身
    delete hash;
}


//
//void insert_pos_exist(bit64Hash* p_root, uint32_t arrayID, uint64_t seq_current, uint32_t pos) {
//	uint32_t shift = bit64HashFunction(seq_current);
//	p_root->p_pos[shift][arrayID].push_back(pos);
//}
