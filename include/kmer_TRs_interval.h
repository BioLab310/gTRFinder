

#ifndef INCLUDE_KMER_TRS_INTERVAL_H_
#define INCLUDE_KMER_TRS_INTERVAL_H_

#include "basic.h"
typedef vector<pair<int,int>> startEnd;
struct kmer_TRs_interval{
	char** kmer;
	vector<startEnd> start_ends;
	vector<int> spaces;
	int length=0;

    void swapElement(int i, int j);
    int partition(int low, int high);
    void quickSort(int low, int high);
    void sortBySpace();
};


#endif /* INCLUDE_KMER_TRS_INTERVAL_H_ */
