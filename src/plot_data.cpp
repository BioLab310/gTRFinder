
#include "plot_data.h"

std::vector<std::pair<uint32_t,uint32_t>> find_boundary_indices(const std::vector<uint32_t>& values, const std::vector<std::pair<int,int>>& boundaries) {
    vector<pair<uint32_t,uint32_t>> index_ranges; // 存储符合条件的索引范围
    for (const auto& boundary : boundaries) {
        int low = boundary.first;   // 边界的下限
        int high = boundary.second;  // 边界的上限
        int start_index = -1;       // 用于记录值的起始索引
        int end_index = -1;         // 用于记录值的结束索引

        for (int i = 0; i < values.size()-2; ++i) {
            if (values[i] >= low && values[i] <= high) {
                if (start_index == -1) {
                    start_index = i; // 记录起始索引
                }
                end_index = i; // 不断更新结束索引
            }
        }
        if (start_index != -1 && end_index != -1) {
            index_ranges.emplace_back(start_index, end_index);
        }
    }
    return index_ranges; // 返回所有找到的索引范围
}
pair<vector<uint32_t>, vector<uint32_t>> true_false_compute(vector<uint32_t> index){
	vector<uint32_t> solid_index;
	vector<uint32_t> hollow_index;
	vector<pair<uint32_t,uint32_t>> boundaries;
	for(int i=0;i<index.size()-2;i++){
		int is_in_bound=false;
		for(auto j:seq_bound){
			if(index[i]<=j.second&&index[i]>j.first){
				is_in_bound=true;
			}
		}
		if(is_in_bound){
			solid_index.push_back(i);
		}
		else{
			hollow_index.push_back(i);
		}
	}
	if(index.size()<2) false_charactor_ratio.push_back(0);
	else{
		false_charactor_ratio.push_back(hollow_index.size()*1.0/(index.size()-2));
	}
	return {solid_index, hollow_index};
}
void plot_data_by_py(vector<double> data,int kmer,int top, vector<uint32_t> index, int flag){
	auto[solid_index, hollow_index]=true_false_compute(index);
	vector<pair<uint32_t,uint32_t>> boundaries;
	boundaries=find_boundary_indices(index,seq_bound);
	if(flag==2) {
		plot_histogram(data, kmer);
	}
	else{
		plot_scatter(data, kmer, top,solid_index ,hollow_index, boundaries, flag);//绘制散点图
	}
//	plot_histogram(data, kmer);
//	cout<<"in generate"<<' '<<flag<<endl;
	return;
}
void plot_scatter(vector<double> data, int kmer,int top,vector<uint32_t> solid_index, vector<uint32_t> hollow_index,
		vector<pair<uint32_t,uint32_t>> boundaries, bool flag){
	try {
	    // 导入 Python 文件
	    py::module sys = py::module::import("sys");
	    sys.attr("path").attr("append")("../src/");//plot_data文件在生成文件前一个文件
//	    py::print(sys.attr("path"));
	    // 准备数据
	    py::module plot_data_module = py::module::import("plot_data");
	    string title = to_string(kmer); // C-style string
	    // 将 C++ vector<int> 转换为 Python 列表
	    py::list py_data;
	    for (const auto& val : data) {
	        py_data.append(val);
	    }
	    py::list py_boundaries;
		for (const auto& boundary : boundaries) {
			py_boundaries.append(py::make_tuple(boundary.first, boundary.second));
		}
		py::list py_solid_index;
		for (const auto& index : solid_index) {
			py_solid_index.append(index);
		}
		py::list py_hollow_index;
		for (const auto& index : hollow_index) {
			py_hollow_index.append(index);
		}
	    // 调用 Python 函数
	    plot_data_module.attr("plot_data")(py_data, title,top,py_solid_index,py_hollow_index,py_boundaries,init_seq_time,
	    		edit_distance_radio_global, flag);

	}
	catch (const py::error_already_set& e) {
	    std::cerr << "Python error: " << e.what() << std::endl;
	}
	catch (const std::exception& e) {
	    std::cerr << "Standard error: " << e.what() << std::endl;
	}
}

void plot_histogram(vector<double> data, int kmer){
	try {
	    // 导入 Python 文件
	    py::module sys = py::module::import("sys");
	    // sys.attr("path").attr("append")("../");//plot_data文件在生成文件前一个文件
	    // 准备数据
	    py::module plot_data_module = py::module::import("plot_data");
	    string title = to_string(kmer); // C-style string
	    // 将 C++ vector<int> 转换为 Python 列表
	    py::list py_data;
	    for (const auto& val : data) {
	        py_data.append(val);
	    }
	    // 调用 Python 函数
	    plot_data_module.attr("plot_histogram")(py_data, init_seq_time,len_global, time_global,edit_distance_radio_global);

	}
	catch (const py::error_already_set& e) {
	    std::cerr << "Python error: " << e.what() << std::endl;
	}
	catch (const std::exception& e) {
	    std::cerr << "Standard error: " << e.what() << std::endl;
	}
}

