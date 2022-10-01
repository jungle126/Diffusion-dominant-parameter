# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 16:44:07 2022

@author: Jungle
"""
import math
import matplotlib.pyplot as plt
from scipy import interpolate

pi = 3.14
def Rs(t,mi):
    a = 4*t
    b = math.log(4*math.sqrt(pi**3*t**3)/(10**mi))
    Rs = math.sqrt(-a*b)
    return Rs
def Rs_t_mi(t_min,t_max,timestep,mi):
    
    listt = [i*timestep for i in range(t_min,t_max+1)]
    listR  = []
    for i in range(len(listt)):
        t = listt[i]
        Rsi = Rs(t,mi)
        listR.append(Rsi)
        print('t = {}s,Rs = {}m'.format(t,Rs))
    return listt,listR

td_list = [1187.1, 1895.5, 4112.2, 6512.7, 9096.0, 11180.4, 12811.2, 14964.0, 17848.8]
Rt_list = [57858.34161380612,100005.1550379589,149908.5177655292,
 220829.67046665933,473711.840161347,752612.6559436886,829971.0862877094,
 869325.0105844229,894994.7187997375]
FigurePath = 'D:\gycx\pic\\'
FigureName = 'Rd-td.jpg'
#plot
fig=plt.figure(figsize=(8,6),dpi=300)#添加画布
x = td_list
y = Rt_list
xnew = [i for i in range(int(x[0])+1,int(x[-1]))]

# plt.plot(x,y,color = 'k',linewidth = 1.2,marker = '^',markersize = 1.2)

# 插值
kind = ["cubic"]
# pltcolor = ['r','g','b','k','p']
for i in range(len(kind)):#插值方式
    #"nearest","zero"为阶梯插值
    #slinear 线性插值
    #"quadratic","cubic" 为2阶、3阶B样条曲线插值
    f=interpolate.interp1d(x,y,kind=kind[i])
    # ‘slinear’, ‘quadratic’ and ‘cubic’ refer to a spline interpolation of first, second or third order)
    ynew=f(xnew)
    #plt.plot(xnew,ynew,linewidth = 1.2,label='Rt')
ynew2 = [] #Rd
for i in xnew: #i是min
    Rdi = Rs(i,10)
    ynew2.append(Rdi)
plt.plot(xnew,ynew2,color = 'r',linewidth = 1.2,label='Rd')
plt.xlim(0, 19600)
plt.xticks([i*1800 for i in range(12)],[i*30 for i in range(12)])
#plt.ylim(0,950000)
plt.xlabel('td(min)',size = 13)
plt.ylabel('Rd(m)', size = 13)
plt.text(500, 900000, 'emit hight = 3000m,drop time = 600min, different sedimentation velocity',fontsize = 8) #text前面的标是以图的刻度为参考的
# plt.text(500, 615, 'emit hight = 3000m,drop time = 600min',fontsize = 8)
plt.legend(loc="lower right")
plt.savefig(FigurePath + FigureName)