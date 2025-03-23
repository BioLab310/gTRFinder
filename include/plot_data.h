

#ifndef INCLUDE_PLOT_DATA_H_
#define INCLUDE_PLOT_DATA_H_
#include "basic.h"
#include <pybind11/embed.h>
namespace py = pybind11;
void plot_data_by_py(vector<double> data, int kmer,int top,vector<uint32_t> index, int flag);
void plot_scatter(vector<double> data, int kmer,int top,vector<uint32_t> solid_index, vector<uint32_t> hollow_index,
		vector<pair<uint32_t,uint32_t>> boundaries, bool flag);
void plot_histogram(vector<double> data, int kmer);
pair<vector<uint32_t>, vector<uint32_t>> true_false_compute(vector<uint32_t> index);


#endif /* INCLUDE_PLOT_DATA_H_ */
