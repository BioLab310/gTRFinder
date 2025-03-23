
#include "evaluate.h"

// 越大越好 - 标准化函数
std::vector<double> normalizeMaximize(const std::vector<double>& X) {
    double X_min = *std::min_element(X.begin(), X.end());
    double X_max = *std::max_element(X.begin(), X.end());

    std::vector<double> X_normalized(X.size());

    for (size_t i = 0; i < X.size(); ++i) {
        if (X_max == X_min) {
            X_normalized[i] = 0.5; // 防止除以零
        } else {
            X_normalized[i] = (X[i] - X_min) / (X_max - X_min);
        }
    }

    return X_normalized;
}

// 越小越好 - 标准化函数
vector<double> normalizeMinimize(const std::vector<double>& Y) {
    double Y_min = *std::min_element(Y.begin(), Y.end());
    double Y_max = *std::max_element(Y.begin(), Y.end());

    std::vector<double> Y_normalized(Y.size());

    for (size_t i = 0; i < Y.size(); ++i) {
        if (Y_max == Y_min) {
            Y_normalized[i] = 0.5; // 防止除以零
        } else {
            Y_normalized[i] = (Y_max - Y[i]) / (Y_max - Y_min);
        }
    }

    return Y_normalized;
}

// 综合评价函数，返回分数和最高分数的索引
uint32_t evaluateScore(const std::vector<double>& X, const std::vector<double>& Y, double w_X, double w_Y) {
    vector<double> X_normalized = normalizeMaximize(X);
    vector<double> Y_normalized = normalizeMinimize(Y);

    std::vector<double> scores(X.size());
    size_t max_index = 0;
    double max_score = -1.0;

    for (uint32_t i = 0; i < X.size(); ++i) {
        scores[i] = w_X * X_normalized[i] + w_Y * Y_normalized[i];
        // 更新最高分和索引
        if (scores[i] > max_score) {
            max_score = scores[i];
            max_index = i;
        }
    }
//    for(int i=0;i<scores.size();i++){
//    	cout<<fixed<<setprecision(3)<<X[i]<<' '<<Y[i]<<' '<<scores[i]<<endl;
//    }

    return max_index+3;  // 返回分数和最高分索引
}

