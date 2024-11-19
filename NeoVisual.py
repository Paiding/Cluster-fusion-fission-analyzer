import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter
import numpy as np

def format_of_plot():
    #plt.xlim(0,10)
    #plt.ylim(0,0.5)
    bwith = 1.5 #边框宽度设置为2
    ax = plt.gca()#获取边框
    #ax.spines['top'].set_color('red')  # 设置上‘脊梁’为红色
    #ax.spines['right'].set_color('none')  # 设置上‘脊梁’为无色
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    
    plt.tick_params(axis='both',which='major',direction='in',labelsize=15,pad=12,length=8,width=1)
    plt.tick_params(axis='both',which='minor',direction='in',labelsize=15,pad=12,length=4,width=1)

    plt.rcParams['font.sans-serif'] = ['Times New Roman']  # 设置字体为SimHei显示中文
    plt.rcParams['axes.unicode_minus'] = False  # 设置正常显示符号
def custom_ticks(x, pos):
    """自定义刻度格式化函数"""
    if x > 100:
        return f'{int((x - 99) * 100)}'
    return int(x)

#########################################读数据###############################################

# 读取文件
text = 'cg_alpha_mix'
with open('split_'+text+'.dat', 'r') as f:
    lines = f.readlines()

# 计算矩阵的大小
n = int(np.sqrt(len(lines)))

# 初始化矩阵
matrix = np.zeros((n, n))

# 填充矩阵
for i in range(n):
    for j in range(n):
        x = float(lines[i*n + j].strip())
        matrix[i][j] = x
        #if(x > 0):
            #matrix[i][j] = np.log10(x)
        #else:
            #matrix[i][j] = -4
dx,dy = 1,1
y, x = np.mgrid[slice(1, 187 + dy, dy),
                slice(1, 187 + dx, dx)]


# 创建一个FuncFormatter对象，使用自定义的刻度格式化函数
formatter = FuncFormatter(custom_ticks)
fig = plt.figure(figsize=(8.0,6.4),dpi=150)
ax = plt.subplot()
cmap = plt.get_cmap('Spectral')
matrix = matrix[:-1, :-1]
levels = MaxNLocator(nbins=15).tick_values(matrix.min(), matrix.max())
# contours are *point* based plots, so convert our bound into point
# centers
#cf = ax1.contourf(x[:-1, :-1] + dx/2.,y[:-1, :-1] + dy/2., matrix, levels=levels,cmap=cmap)
#cf = ax1.contourf(x[:-1, :-1] + dx/2.,y[:-1, :-1] + dy/2., matrix, norm=LogNorm(vmin=matrix.min(),vmax=matrix.max()),cmap=cmap)
#cf = ax1.contourf(x[:-1, :-1] + dx/2.,y[:-1, :-1] + dy/2., matrix, norm=LogNorm(),cmap=cmap)
cf = ax.pcolor(x[:-1, :-1] + dx/2.,y[:-1, :-1] + dy/2., matrix, norm=LogNorm(),cmap=cmap,shading='flat')
cbar = fig.colorbar(cf, ax=ax)
#设置colorbar上的字体大小
cbar.ax.tick_params(labelsize=14)
#ax1.set_title('contourf with levels')

format_of_plot()

font_size = {'size':28}
ax.xaxis.set_major_formatter(formatter)
ax.yaxis.set_major_formatter(formatter)
#plt.title('Rate-Cluster size',font_size,pad=20)
plt.xlabel('Cluster Aggregation Number',fontsize=16)
plt.ylabel('Cluster Aggregation Number',fontsize=16)
#plt.legend(loc='best',fontsize=20)

# adjust spacing between subplots so `ax1` title and `ax0` tick labels
# don't overlap
fig.tight_layout()
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.savefig(text + '.png')
plt.show()
