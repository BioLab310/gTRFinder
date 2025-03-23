
#ifndef HASH_H_
#define HASH_H_
#include "basic.h"
#define HashSize 0x100000
#define MinSizeForSearch 100
#define AddSize 1
struct bit64Hash {
	uint64_t** p_kmer;//k_mer的值,hash原值
	vector<vector<vector<uint32_t>>> p_pos;//存放k_mer位置
	uint64_t* len;
	uint64_t* maxlen;
};


struct bit64Hash* initial256BitHashTable();
uint32_t bit64HashFunction(uint64_t current);
int32_t search64BitHashTable(bit64Hash* ph, uint64_t seq_current);
int32_t binary_64_bit_search(vector<vector<uint32_t>> pos, uint64_t* p, uint64_t len, uint64_t key);
void insert64BitHashTable(bit64Hash* ph, uint64_t kmer, uint32_t pos);
void insert_pos_exist(bit64Hash* p_root, uint32_t arrayID, uint64_t seq_current, uint32_t pos);
void deleteBit64Hash(bit64Hash* hash, uint64_t kmer_size);




#endif /* HASH_H_ */
