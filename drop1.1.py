import netCDF4 as nc
import numpy as np
import bisect
import math
import matplotlib.pyplot as plt
import re
import os 
import matplotlib.gridspec as gridspec
import seaborn as sns

# 距离转换成经纬度插值算法（x，y）单位m->（经度，纬度）
def xy_latlong(x,y,sstep):
    '''
    x:单位m
    y:单位m
    sstep：距离步长，wrf数据中一格为多少m，
    '''
    x_238_index = x/sstep
    y_188_index = y/sstep
    x_238_index_min = int(x_238_index)
    y_188_index_min = int(y_188_index)
    
    w00 = abs((1-(x_238_index-x_238_index))*(1-(y_188_index-y_188_index_min)))
    w01 = abs((1-(x_238_index-x_238_index))*(y_188_index-y_188_index_min))
    w10 = abs((x_238_index-x_238_index)*(1-(y_188_index-y_188_index_min)))
    w11 = abs((x_238_index-x_238_index)*(y_188_index-y_188_index_min))
    
    LAT00 = df.variables['XLAT'][0][x_238_index_min][y_188_index_min]
    LAT01 = df.variables['XLAT'][0][x_238_index_min][y_188_index_min + 1]
    LAT10 = df.variables['XLAT'][0][x_238_index_min + 1][y_188_index_min]
    LAT11 = df.variables['XLAT'][0][x_238_index_min + 1][y_188_index_min + 1]
    
    LONG00 = df.variables['XLONG'][0][x_238_index_min][y_188_index_min]
    LONG01 = df.variables['XLONG'][0][x_238_index_min][y_188_index_min + 1]
    LONG10 = df.variables['XLONG'][0][x_238_index_min + 1][y_188_index_min]
    LONG11 = df.variables['XLONG'][0][x_238_index_min + 1][y_188_index_min + 1]
    
    XLAT = w00*LAT00+w01*LAT01+w10*LAT10+w11*LAT11
    YLONG = w00*LONG00+w01*LONG01+w10*LONG10+w11*LONG11
    return XLAT, YLONG
    
    

#最外层嵌套



#### function
# 距离转经纬度

# 落点边界确定，包含所有落点的最小矩形，向四周扩5000m，每100m一个计算点
def set_boundry(listx,listy):
    '''输入所有点的xy坐标，所有x坐标放在listx中，所有y坐标放在listy中，不需要一一对应
        返回边界list和四个角（单位是新建立的格点，一格100m）    
    '''
    x_min = min(listx)*sstep//100*100-6000
    x_max = max(listx)*sstep//100*100+6000
    y_min = min(listy)*sstep//100*100-6000
    y_max = max(listy)*sstep//100*100+6000
    x_sep = int((x_max-x_min)//100+1)
    y_sep = int((y_max-y_min)//100+1)
    boundary = [[0 for i in range(x_sep)] for j in range(y_sep)]
    return x_min,x_max,y_min,y_max,boundary

# 框定单个落点所需计算的周围 5000m内的范围，101*101个点
def coculate_point_suround(X,Y,x_min,y_min):
    index_X,index_Y = XY_index(X,Y,x_min,y_min)
    index_Xmin = index_X-50
    index_Xmax = index_X+50
    index_Ymin = index_Y-50
    index_Ymax = index_Y+50
    return index_Xmin,index_Xmax,index_Ymin,index_Ymax

#    for i in range(index_Ymin,index_Ymax+1):
#        for j in range(index_Xmin,index_Xmax+1):
#            X0,Y0 = gridindex_XY(j,i,x_min,y_min)
#            r = R_func(X0, Y0, H0, Xs, Ys, Zs)
#            boudary[i][j] += Q_sim(r,sigma)

# 格点的index_y和index_x转换成Y与X    
def gridindex_XY(index_X,index_Y,x_min,y_min):
    X = x_min+100*index_X
    Y = y_min+100*index_Y
    return X,Y

# XY转格点
def XY_index(X,Y,x_min,y_min):
    index_X = int((X-x_min)//100-1)
    index_Y = int((Y-y_min)//100-1)
    return index_X,index_Y
    

# 浓度计算简化版，只要知道距离排放源落地点与所求点的距离R以及sigma的值来得到结果
def Q_sim(r,sigma,pi = 3.14):
    Q = math.exp(-r**2/(2*sigma**2))/(math.sqrt(2)*pi**3*sigma**3)
    return Q

# 一个用来计算R的func,仅用于计算落点和所求点距离
def R_func(X0,Y0,Z0,Xs,Ys,Zs):#所求点的X0Y0Z0，落点的XsYsZs，单位的的都为m
    r_xy = math.sqrt((Xs-X0)**2+(Ys-Y0)**2)
    delt_z = abs(Z0-Zs)
    r = math.sqrt(r_xy**2+delt_z**2)
    return r

# 一个用来计算sigma的func
def sigma_func(t,D = 1): #t是落地所经历的时间，D是扩散参数，这里默认取1
    sigma = math.sqrt(2*D*t)
    return sigma

# 用来对浓度进行对数处理
def log_10(value):
    if value:
        k = math.log10(value)
        if k > -200:
            return k
        else:
            return -200
    else:
        return -200
##########
sstep = 27000  #网格步长
df = nc.Dataset("D:/gycx/zzytest.nc")
df_U = np.array(df['U'])
df_V = np.array(df['V'])
df_W = np.array(df['W'])
df_HGT = np.array(df["HGT"])
df_H = np.array(df["H"])
lists = []
set_height = [3000]
#set_height = [3000,5000,8000]
set_V = [1.1]
set_V = [0.249,0.26,0.27,0.28,0.3,0.35,0.4,0.5,0.7,1.1,2.2,3.5]
set_t = np.arange(0,60,1) # （a,b,c）a表示第一颗从什么时候开始排（单位分钟），b排到啥时候，c间隔
id_list = list(range(len(set_t))) # 用于屏幕输出第几次排放
set_m = ['Jan','Apr','Jul','Oct']
for sh in set_height: 
    for sv in set_V:
        x_y_z_t_sigma = [[] for i in range(5)] # sh，sv排列组合，生成多张图
        for st in range(len(set_t)):
            
            num = 0    ##单纯记录次数
            T = 3600   ##数据分辨率小时
            ## 时间记录为t_h:t_min
            t_h = set_t[st]//60            ##实时小时
            t_min = set_t[st] - t_h*60     ##实时时间 单位分钟
            
            x,x_title = 80,80     ##初始点坐标 这里是格点，238*188
            y,y_title = 70,70
            #z = 9445.736
            z = sh + df_HGT[1][int(y)][int(x)]  #排放源的高度，单位为m
            #x,y,z是排放源初始点 ，xy都是格点，z是m
            g = sv    ## 下降速度
            tstep = 60  ##时间步长，单位s
            print('{}：{}时刻 距地高度处{}m，下落速度为{}，第{}个排放开始'.format(t_h,t_min,sh,sv,id_list[st]+1))
            print()

            tt=[]  ##
            ww=[]
            xx=[]
            yy=[]
            suma=[]
            sumb=[]
            #第一次进入时生成H高度数组，每过一个小时更新一次H
            h0 = 0
            H = []
            for i in range(df_H.shape[1]):
                z0 = df_H[t_h][i][int(y)][int(x)]
                h0 = h0 + np.array(z0)  #累加高度h0
                H = np.append(H,h0)
            H = H.tolist()
            print('在开始时形成一个高度 H list')
            print()
            while z > df_HGT[1][int(y)][int(x)]:
                #print('{}:{}时刻排放源中心高度为{}m'.format(t_h,t_min,abs(z-0)))
                if t_min//60: #过一个小时t_min归零，t_h加一
                    print()
                    print("形成一个新的 H list！此后一小时的高度层用新数据")
                    print()
                    t_min = 0
                    t_h += 1
                    h0 = 0
                    H = []  ### 高度集合
                    ### 生成H高度数组
                    for i in range(len(df_H[0])):
                        z0 = df_H[t_h][i][int(y)][int(x)]
                        h0 = h0 + np.array(z0)  #累加高度h0
                        H = np.append(H,h0)
                    H = H.tolist()
                    #print('H数组：',H)
                    #print(abs(z - 0) )
                    
                    ###查询z高度在H高度的层数:简化，向上取整，如3和4层之间取第四层的U V W
                position = bisect.bisect(H, z) 
                #print(t,position,x,y) # 如果放到有序序列中，应该存在的索引位置
                #print("h:",t_real,"s:",num,"总时间",num+t_real*3600,"层数：",position,"x(km):",x,"y(km):",y) # 如果放到有序序列中，应该存在的索引位置

                x_0 = math.floor(x)
                y_0 = math.floor(y)

                xyz_0 = [x_0,x_0+1]
                xyz_1 = [y_0,y_0+1]
                xyz_2 = [position,position+1]
                #print(xyz_0,xyz_1,xyz_2)
                U = []
                V = []
                W = []
                
                #遍历八个格点权重插值
                #print(str(num+t_real*3600)+"秒时的坐标点为")
                for p in xyz_0:
                    for q in xyz_1:
                        for k in xyz_2:
                            #print('坐标点：',p,q,k)
                            b = abs((p-x)*(q-y)*(z-(H[k]))/(H[position+1]-H[position]))
                            U.append(math.sqrt((df_U[t_h][k][q][p])**2*b)/abs(df_U[t_h][k][q][p])*df_U[t_h][k][q][p])
                            V.append(math.sqrt((df_V[t_h][k][q][p])**2*b)/abs(df_V[t_h][k][q][p])*df_V[t_h][k][q][p])
                            W.append(df_W[t_h][k][q][p]*b)
                u = sum(U)
                v = sum(V)
                w = sum(W)
                ws=np.sqrt(u*u+v*v)
                
                #print(str(num+t_real*3600)+"秒时的风速为",(u,v,w))  
                #print(ws)
                x1=[]
                y1=[]
                ###模拟运动后下一秒的位置结果
                # t = t+1
                x = x + v*0.001/3*tstep
                y = y + u*0.001/3*tstep
                z = z + (w - g)*tstep
                
                # if (num+t_real*3600)%3600 == 0:
                #print("当t=",num + t_real*3600 ,"此时坐标：",x,y,z)
                #print("\n")
                t_min += 1
                #num = num + tstep 
                
                tt=np.append(tt,t_h*3600+t_min*60)
                ww=np.append(ww,ws)
                xx=np.append(xx,x)
                yy=np.append(yy,y)
            #确认落点x0,y0，落地时刻t0,落点高度z0
            x0 = xx[len(xx)-1]  # 单位 格点
            y0 = yy[len(yy)-1]  # 单位 格点
            t0 = tt[len(tt)-1]  # 单位 s秒 
            z0 = df_HGT[1][int(y)][int(x)] # 单位 m
            ts = tt[-1] - tt[0] # 落地历时 单位 s
            sigma = sigma_func(ts)
            #存储每一个落点的x，y，z，xy是格点，z是m；存储所经时间（s），计算所得的sigma
            x_y_z_t_sigma[0].append(x0) 
            x_y_z_t_sigma[1].append(y0)
            x_y_z_t_sigma[2].append(z0)
            #x_y_z_t_sigma[3].append(H)
            x_y_z_t_sigma[3].append(ts)
            x_y_z_t_sigma[4].append(sigma)
       
       