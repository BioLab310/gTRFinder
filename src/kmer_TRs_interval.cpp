

#include "kmer_TRs_interval.h"


void kmer_TRs_interval::swapElement(int i, int j){
	 if (i < 0 || j < 0 || i >= length || j >= length) {
		throw out_of_range("Index out of range");
	}

	// swap(kmer[i], kmer[j]);
	swap(start_ends[i], start_ends[j]);
	swap(spaces[i], spaces[j]);
}
int kmer_TRs_interval::partition(int low, int high){
	int pivot=spaces[high];
	int i=low-1;
	for(int j=low; j<high; j++){
		if(spaces[j]< pivot){
			i++;
			swapElement(i, j);
		}
	}
	swapElement(i+1, high);
	return i+1;
}

void kmer_TRs_interval::quickSort(int low, int high){
	if(low<high){
		int p=partition(low, high);
		quickSort(low, p-1);
		quickSort(p+1, high);
	}
}

void kmer_TRs_interval::sortBySpace(){
	if(length>1){
		quickSort(0, length-1);
	}
}



