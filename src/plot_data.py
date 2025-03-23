import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.ticker import LogLocator
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import MultipleLocator
from fractions import Fraction
from matplotlib import font_manager
from matplotlib.ticker import AutoLocator,AutoMinorLocator
def plot_data(y, kmer,top, solid_index, hollow_index,boundaries,init_seq_times, edit_distance_radio, flag):
    if flag==1:
        flag='best'
    else:
        flag=''
    x = np.arange(len(y))  
    y = np.array(y) 
    y_fold = np.where(y < 0.75, 1 / y, y)  
    colors = ['#1E90FF' if value > 0.75 else '#FF4500' for value in y]
    colors=np.array(colors)
    
    maxlim=10
    y_ticks_above = np.arange(1, maxlim, 1)  # 1 以上的刻度（1到10）
    y_ticks_below = 1 / y_ticks_above  # 1 以下的刻度（1 以上的倒数）
    def to_fraction(x, pos):
        # 使用 Fraction 自动将小数转换为分数
        fraction = Fraction(x).limit_denominator()  # 这里可以指定最大分母
        return str(fraction)  # 返回分数字符串
    # 应用自定义格式化函数
    formatter = FuncFormatter(to_fraction)
    #绘制全局变量
    plt.figure(figsize=(7, 7))
    plt.rcParams['font.size'] = 20  # 默认字体大小
    plt.rcParams['font.family']='Times New Roman'
    
    # 绘制原点
    plt.scatter(x[solid_index], y_fold[solid_index], alpha=0.6, color=colors[solid_index], edgecolors='w', linewidths=0.5, marker='o', label='right radio')
    plt.scatter(x[hollow_index], y_fold[hollow_index], alpha=0.6, facecolors='none', edgecolors='red', linewidths=0.5, marker='o', label='false radio')
    
   # 绘制边界
    for boundary in boundaries:
        plt.axvline(x=boundary[0], color='gray', linestyle='--', linewidth=1)
        plt.axvline(x=boundary[1], color='gray', linestyle='--', linewidth=1)
    plt.gca().yaxis.set_major_formatter(formatter)
    #plt.title(f'Pattern Feature')
    
    #绘制x,y轴
    plt.ylim(0,maxlim)
    y_ticks =  y_ticks_above
    plt.yticks(y_ticks, labels=y_ticks)
    plt.xlabel('Position',labelpad=8)
    plt.ylabel('Pattern Feature',labelpad=10)
    plt.tick_params(axis='x', pad=8)  # X 轴数字与坐标轴的距离
    plt.tick_params(axis='y', pad=15)  # Y 轴数字与坐标轴的距离
    
    #设置网格
    plt.grid(True)
    plt.gca().xaxis.set_major_locator(AutoLocator())  # 主刻度每隔1
    # 设置次刻度间隔
    plt.gca().xaxis.set_minor_locator(AutoMinorLocator())  # 次刻度每隔0.2
    plt.gca().set_autoscalex_on(True)  # 自动调整横坐标
    plt.gca().set_autoscaley_on(False)  # 禁止纵坐标自动调整
    plt.grid(which='major', linestyle='-', alpha=0.6)
    plt.grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.6, color='gray')
    
    plt.gca().spines['top'].set_alpha(0.4)     # 上边框
    plt.gca().spines['right'].set_alpha(0.4)   # 右边框
    plt.gca().spines['bottom'].set_alpha(0.4)  # 下边框
    plt.gca().spines['left'].set_alpha(0.4)    # 左边框
    #plt.legend()  # 显示图例
    plt.savefig(fr'scatter/{init_seq_times}_{kmer}{flag}_{edit_distance_radio:.2f}_ratio_scatter.jpg')  # 保存图像
    plt.close()
def custom_mapping(x):
    # 映射规则
    if x == 0.8:
        return -0.5  # 将 0.5-1 的区间拉伸
    elif x==1.25:
        return 0.5;
    if x<1:
        return -1/x+0.5 
    else :
        return x-0.5  
def plot_histogram(data, init_seq_time,lens, time,edit_distance_radio):
    data = np.array(data) 
    # 设置左右边界的间隔
    n= 1
    bins_right = np.array([1.25,2,3,4,5,6,7,8,9,10])
    bins_left=sorted(1/bins_right)
    bins=np.concatenate((bins_left, bins_right))
    bin_map=np.empty(bins.size)
    for i in range(bins.size):
        bin_map[i]=custom_mapping(bins[i])
    # 绘制直方图
    bin_width = np.diff(bins).min()  # Assuming all bins have a similar width
    counts,_=np.histogram(data, bins=bins)
    # 设置自定义的刻度和标签
    custom_ticks = bins  # 自定义刻度
    custom_labels = [f"{tick:.3g}" for tick in bins]  # 将界限转换为文本
    def to_fraction(x):
        # 使用 Fraction 将小数转换为分数
        fraction = Fraction(x).limit_denominator()  # 自动转为最简分数形式
        numerator, denominator = fraction.numerator, fraction.denominator
        if denominator == 1:  # 如果是整数，直接显示
            return numerator
        else:  # 否则显示分数形式
            return rf'$\frac{{{numerator}}}{{{denominator}}}$'
    font_path = '/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf'  # 确保路径正确
    prop = font_manager.FontProperties(fname=font_path)
    
    plt.figure(figsize=(7, 7))
    plt.rcParams['font.size'] = 15 # 默认字体大小
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman']
    plt.bar(np.array(bin_map[:-1])+0.5, counts, width=1, 
        color='skyblue', edgecolor='black', linewidth=0.8, alpha=1)
    
 #   plt.gca().xaxis.set_major_locator(LogLocator(base=1.2)) 
    #plt.title(f'Distribution of Radio Values for TR Features of {init_seq_time}_{lens}_{time}_{edit_distance_radio:.2f}')
    plt.xlabel('Ratio Values', labelpad=6,fontsize=20)
    plt.ylabel('Frequency', labelpad=8, fontsize=20)
    plt.xticks(bin_map,labels= np.vectorize(to_fraction)(bins), fontproperties=prop)
    plt.grid(True, which="both", ls="--", lw=0.5,  alpha=0.6, color='gray')
    plt.axvline(x=0, color='#FF4500', linestyle='--', linewidth=1,label='x = 1')
# 在虚线旁边添加标签
    #plt.text(0.25, max(counts)*1.01, 'x=1', color='red', fontsize=12, horizontalalignment='center')
    plt.legend(fontsize=15)
    plt.savefig(fr'histogram/{init_seq_time}_{lens}_{time}_{edit_distance_radio:.2f}_ratio_histogram.svg')  # 保存图像
    plt.close()

