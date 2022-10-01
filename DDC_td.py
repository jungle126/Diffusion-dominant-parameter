# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 16:44:07 2022

@author: Jungle
"""
import math
import matplotlib.pyplot as plt
from scipy import interpolate


def DDC3_fun(Rd,Rt):
    DDC3 = 1-Rt**2/(Rt+Rd)**2
    return DDC3
def DDC_func(Rd,Rt,te=1,delta_te=1):
    above = Rd**2*(1+te//delta_te)
    below = (Rd+Rt)**2
    DDC = above/below
    
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

x = td_list
y = Rt_list

xnew = [i for i in range(int(x[0])+1,int(x[-1]))] #定义x轴插值范围
Rd_listnew = [] #Rd的赋值
for i in xnew: #i是min
    Rdi = Rs(i,10)
    Rd_listnew.append(Rdi)
# plt.plot(x,y,color = 'k',linewidth = 1.2,marker = '^',markersize = 1.2)

# 插值 得到对于x轴的所有Rt，放在Rt_listnew内
kind = ["cubic"]
# pltcolor = ['r','g','b','k','p']
for i in range(len(kind)):#插值方式
    #"nearest","zero"为阶梯插值
    #slinear 线性插值
    #"quadratic","cubic" 为2阶、3阶B样条曲线插值
    f=interpolate.interp1d(x,y,kind=kind[i])
    Rt_listnew=f(xnew)
   
# 赋值DDC1，DDC2，DDC3
DDC1_list = []
DDC2_list = []
DDC3_list = []
for i in range(len(xnew)):
    DDC1 = DDC_func(Rd_listnew[i],Rt_listnew[i],600)
    DDC2 = DDC_func(Rd_listnew[i],Rt_listnew[i])
    DDC3 = DDC3_fun(Rd_listnew[i],Rt_listnew[i])
    DDC1_list.append(DDC1)
    DDC2_list.append(DDC2)
    DDC3_list.append(DDC3)


# #plot
fig=plt.figure(figsize=(8,6),dpi=300)#添加画布
FigurePath = 'D:\gycx\pic\\'
FigureName = 'DDC—td 100%.jpg'
# plt.text(100, 0.96, 'emit hight = 3000m, emit time = 600min,different sedimentation velocity')
x = xnew
y1 = DDC1_list
y2 = DDC2_list
y3 = DDC3_list
plt.plot(x,y1,color = 'g',linewidth = 1.2,label='DDC1')
plt.plot(x,y2,color = 'r',linewidth = 1.2,label='DDC2')
plt.plot(x,y3,color = 'b',linewidth = 1.2,label='DDC3')

plt.xlim(0, 19600)
plt.xticks([i*1800 for i in range(12)],[i*30 for i in range(12)])
# plt.ylim(0,0.01)
# plt.yticks([0,0.0025,0.0050,0.0075,0.01],[0,0.25,0.5,0.75,1])
plt.ylim(0,1)
plt.yticks([0,0.25,0.50,0.75,1],[0,25,50,75,100])
plt.xlabel('td(min)',size = 13)
plt.ylabel('ratio(%)', size = 13)
plt.legend(loc="lower right")

plt.savefig(FigurePath + FigureName)