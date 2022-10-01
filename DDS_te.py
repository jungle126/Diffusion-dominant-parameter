import math
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np

FigurePath = 'D:\gycx\pic\\'
FigureName = 'DDI-te.jpg'
##Rfunc
def R_func(X0,Y0,Z0,Xs,Ys,Zs):#所求点的X0Y0Z0，落点的XsYsZs，单位的的都为m
    r_xy = math.sqrt((Xs-X0)**2+(Ys-Y0)**2)
    delt_z = abs(Z0-Zs)
    r = math.sqrt(r_xy**2+delt_z**2)
    return r
def DDC3_fun(Rd,Rt):
    DDC3 = 1-Rt**2/(Rt+Rd)**2
    return DDC3
def DDC_func(Rd,Rt,te=1,delta_te=1):
    above = Rd**2*(1+te//delta_te)
    below = (Rd+Rt)**2
    DDC = above/below
    return DDC

sstep = 27000

#sh = 3000 sv = 1.1 
x_y_z_t_sigma = [[95.04160770716841, 95.00529789775035, 95.01919370559729, 95.0340465018938, 95.05846989816953, 94.97303158319816, 94.96939465772343, 95.05714010481401, 94.98508028711383, 94.9982403258967, 94.96442365982378, 95.00763557440571, 94.9863619762275, 94.99624856714912, 95.01915471894962, 94.97051889947849, 94.95820764152442, 94.98523714792407, 95.05807724405314, 94.98408831076628, 95.01751101963647, 94.99339617550174, 95.02482730001518, 94.98428568241395, 95.00862373117833, 95.0219592785889, 94.99211921587016, 94.97629527151138, 94.97628788925554, 95.05495247294338, 95.05651490672915, 94.9931178252508, 94.99386027298964, 95.04776032839109, 94.97614588323083, 95.00980354653903, 94.97217591849852, 94.98219041437707, 95.0533233962987, 94.97498014428945, 94.96987727332741, 94.97515373246502, 95.01852336630334, 94.9749648948904, 95.06073399525584, 95.03487222002127, 95.06259217548113, 94.97604518273029, 94.97594918093549, 95.02735732400099, 95.07850065152931, 95.02674560259284, 94.97756419696375, 94.97014009768021, 95.0956779627379, 95.0986888935962, 95.1028436786599, 95.06194747199574, 95.14962177621075, 95.15333312919059], [72.68478950743376, 72.71993857570442, 72.68752801395993, 72.65016812879765, 72.6376590418842, 72.61458841451545, 72.61111405857305, 72.6004788007866, 72.65621400561596, 72.6379427816354, 72.62196990143663, 72.58974699192613, 72.45141522076865, 72.59446724585506, 72.36777600828282, 72.34267203860078, 72.38639446959833, 72.52786956796133, 72.55805109400582, 72.49941142228579, 72.6721698874305, 72.71149520347727, 72.70745619138147, 72.53372961639673, 72.4077160233407, 72.50035351281248, 72.54884097987554, 72.6135510221427, 72.66926766479692, 72.48880442379445, 72.44807111471478, 72.58300885288182, 72.46385592043582, 72.37912284655786, 72.48773107072316, 72.54357021837704, 72.49102881545328, 72.55120641488807, 72.34251706934788, 72.39105295278345, 72.45580471171203, 72.63396988215052, 72.35321830194702, 72.44299160193013, 72.50750999923, 72.23468886081218, 72.6194135349794, 72.53704769976092, 72.66310714499792, 72.6301073360547, 72.54704705323469, 72.66687417688017, 72.82354657941683, 72.69129976102586, 72.66394435235638, 72.64288169777903, 72.63215178809436, 72.81881021954628, 72.84231206192887, 72.68171530849179], [479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 624.5753, 624.5753, 479.29004, 624.5753, 624.5753, 624.5753, 479.29004, 624.5753, 624.5753, 479.29004, 624.5753, 624.5753, 624.5753, 479.29004, 624.5753, 479.29004, 624.5753, 479.29004, 624.5753, 479.29004, 479.29004, 624.5753, 624.5753, 624.5753, 479.29004, 479.29004, 624.5753, 624.5753, 479.29004, 624.5753, 479.29004, 624.5753, 624.5753, 479.29004, 624.5753, 624.5753, 624.5753, 479.29004, 624.5753, 479.29004, 479.29004, 479.29004, 624.5753, 624.5753, 479.29004, 479.29004, 479.29004, 624.5753, 624.5753, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004], [3840.0, 3840.0, 3840.0, 3840.0, 3840.0, 3720.0, 3720.0, 3840.0, 3840.0, 3840.0, 3720.0, 3840.0, 3840.0, 3780.0, 3840.0, 3780.0, 3840.0, 3840.0, 3840.0, 3720.0, 3840.0, 3780.0, 3840.0, 3780.0, 3840.0, 3840.0, 3780.0, 3720.0, 3840.0, 3840.0, 3840.0, 3780.0, 3720.0, 3840.0, 3720.0, 3840.0, 3720.0, 3720.0, 3840.0, 3720.0, 3840.0, 3840.0, 3840.0, 3780.0, 3840.0, 3840.0, 3840.0, 3780.0, 3840.0, 3840.0, 3840.0, 3840.0, 3840.0, 3780.0, 3840.0, 3840.0, 3840.0, 3840.0, 3840.0, 3840.0], [87.63560920082658, 87.63560920082658, 87.63560920082658, 87.63560920082658, 87.63560920082658, 86.2554346113913, 86.2554346113913, 87.63560920082658, 87.63560920082658, 87.63560920082658, 86.2554346113913, 87.63560920082658, 87.63560920082658, 86.94826047713663, 87.63560920082658, 86.94826047713663, 87.63560920082658, 87.63560920082658, 87.63560920082658, 86.2554346113913, 87.63560920082658, 86.94826047713663, 87.63560920082658, 86.94826047713663, 87.63560920082658, 87.63560920082658, 86.94826047713663, 86.2554346113913, 87.63560920082658, 87.63560920082658, 87.63560920082658, 86.94826047713663, 86.2554346113913, 87.63560920082658, 86.2554346113913, 87.63560920082658, 86.2554346113913, 86.2554346113913, 87.63560920082658, 86.2554346113913, 87.63560920082658, 87.63560920082658, 87.63560920082658, 86.94826047713663, 87.63560920082658, 87.63560920082658, 87.63560920082658, 86.94826047713663, 87.63560920082658, 87.63560920082658, 87.63560920082658, 87.63560920082658, 87.63560920082658, 86.94826047713663, 87.63560920082658, 87.63560920082658, 87.63560920082658, 87.63560920082658, 87.63560920082658, 87.63560920082658]]

#sh = 3000 sv=  1.1 180min
#x_y_z_t_sigma = [[95.04160770716841, 95.00529789775035, 95.01919370559729, 95.0340465018938, 95.05846989816953, 94.97303158319816, 94.96939465772343, 95.05714010481401, 94.98508028711383, 94.9982403258967, 94.96442365982378, 95.00763557440571, 94.9863619762275, 94.99624856714912, 95.01915471894962, 94.97051889947849, 94.95820764152442, 94.98523714792407, 95.05807724405314, 94.98408831076628, 95.01751101963647, 94.99339617550174, 95.02482730001518, 94.98428568241395, 95.00862373117833, 95.0219592785889, 94.99211921587016, 94.97629527151138, 94.97628788925554, 95.05495247294338, 95.05651490672915, 94.9931178252508, 94.99386027298964, 95.04776032839109, 94.97614588323083, 95.00980354653903, 94.97217591849852, 94.98219041437707, 95.0533233962987, 94.97498014428945, 94.96987727332741, 94.97515373246502, 95.01852336630334, 94.9749648948904, 95.06073399525584, 95.03487222002127, 95.06259217548113, 94.97604518273029, 94.97594918093549, 95.02735732400099, 95.07850065152931, 95.02674560259284, 94.97756419696375, 94.97014009768021, 95.0956779627379, 95.0986888935962, 95.1028436786599, 95.06194747199574, 95.14962177621075, 95.15333312919059, 95.08508954426212, 95.13744054924403, 95.18387305213146, 95.225137309432, 95.27019301560979, 95.27840163210206, 95.37390925929257, 95.24603395762917, 95.14529364632729, 95.23074503739564, 95.4163830912856, 95.21176715124906, 95.34954373558462, 95.21444983986882, 95.35989814149269, 95.21678717875783, 95.44518483293484, 95.6983316148452, 95.81039812542, 95.76816526674531, 95.31330898766066, 95.65588490552466, 95.83751774468668, 95.9543601382851, 95.88573680893121, 95.5678050668028, 95.99929327224818, 95.99119825209019, 95.98935170417948, 95.97612260222975, 95.97204025404564, 95.96737175713747, 95.99174870466531, 95.97241754088513, 95.96234723254429, 95.97838074397505, 95.9440563552024, 95.91810797084673, 95.97963724925314, 95.9924580543033, 95.73499086466428, 95.25629166474913, 95.98402046597128, 95.98699820027454, 95.69771094249228, 95.99777531125225, 95.49236615330022, 95.73890246375647, 95.95828771650349, 95.79254596398849, 95.97957306872347, 95.89386073959922, 95.71527479965076, 95.958006444585, 95.99011867918924, 95.98257797714976, 95.9838299434475, 95.98024332018316, 95.99507854981479, 95.99219666466244, 96.0069934014866, 96.017278360962, 95.98880982785201, 96.00496689367216, 95.98024448872297, 95.98991597183472, 95.96930827775547, 95.97536763527428, 95.97331463140216, 96.02490118625009, 95.9642388944642, 95.99009605385686, 96.01297724315657, 96.01751799266253, 96.02151635349601, 95.97474258967138, 95.9787176801672, 95.97586833309269, 95.97224178861366, 95.98106969702802, 95.98172328409561, 95.9785831649353, 95.99893955910152, 95.9736477508896, 95.96928899862368, 95.97163803738694, 95.97146961857564, 95.97044168259664, 95.98436825955956, 95.98548825620877, 95.97005350159432, 95.97154830648375, 95.98333521286554, 95.99074416674382, 95.99755059348175, 95.98308567253864, 95.98847824797687, 95.98236335095709, 96.11690240066109, 96.10527380954936, 96.10226572570366, 96.11080035468683, 95.99311443681607, 95.99749341135487, 96.09046546364499, 96.11393410946026, 95.99429161383712, 96.09575256001553, 96.11192158152807, 96.07327478053728, 95.99782431527044, 95.96527933826367, 96.0329370372221, 95.9910951779576, 96.12671396933607, 96.16091140964379, 96.13163056616207, 96.08862623632956, 96.22119531153099, 96.17801425465768, 96.1567699032178], [72.68478950743376, 72.71993857570442, 72.68752801395993, 72.65016812879765, 72.6376590418842, 72.61458841451545, 72.61111405857305, 72.6004788007866, 72.65621400561596, 72.6379427816354, 72.62196990143663, 72.58974699192613, 72.45141522076865, 72.59446724585506, 72.36777600828282, 72.34267203860078, 72.38639446959833, 72.52786956796133, 72.55805109400582, 72.49941142228579, 72.6721698874305, 72.71149520347727, 72.70745619138147, 72.53372961639673, 72.4077160233407, 72.50035351281248, 72.54884097987554, 72.6135510221427, 72.66926766479692, 72.48880442379445, 72.44807111471478, 72.58300885288182, 72.46385592043582, 72.37912284655786, 72.48773107072316, 72.54357021837704, 72.49102881545328, 72.55120641488807, 72.34251706934788, 72.39105295278345, 72.45580471171203, 72.63396988215052, 72.35321830194702, 72.44299160193013, 72.50750999923, 72.23468886081218, 72.6194135349794, 72.53704769976092, 72.66310714499792, 72.6301073360547, 72.54704705323469, 72.66687417688017, 72.82354657941683, 72.69129976102586, 72.66394435235638, 72.64288169777903, 72.63215178809436, 72.81881021954628, 72.84231206192887, 72.68171530849179, 72.87068486759715, 72.82016075677495, 72.7676148327942, 72.73351603458826, 72.69444919963418, 72.66814861727696, 72.58920544903806, 72.72572657347435, 72.83202036844405, 72.75037619670988, 72.56595563726707, 72.66590375408131, 72.50460808994538, 72.57338372957821, 72.69999042901746, 72.19584548377793, 72.07523696273421, 71.90046911029046, 71.80394738650287, 71.7922512669776, 72.70865203383224, 72.42999366464029, 72.16507058063547, 72.06234371165628, 72.04615886800504, 72.3492481629077, 72.05106478200032, 71.99517456824296, 71.97883164726277, 71.91705421991279, 71.91899632930712, 71.87839058600339, 71.73697195356095, 71.86295329567095, 71.82558111756929, 71.65986435529379, 72.0966924922437, 71.16135922365196, 71.93021079134961, 71.86037992517393, 71.42527038669066, 71.6599039018294, 71.17416238946524, 71.12857296428867, 71.36483675073411, 71.08818602457544, 71.61671193457856, 72.17050290414628, 72.01036149169428, 72.10669974045601, 72.05051585640588, 72.03395696821947, 72.136094439479, 72.05249347912648, 71.98745129823254, 72.21560775259269, 72.25576272361836, 72.2756144287668, 72.28675964340104, 72.31245498825432, 72.27030887031644, 72.25450344744938, 72.24004859308602, 72.26770573793716, 72.25104689634371, 72.28176634238368, 72.22087281545704, 72.27955667555479, 72.25504649192773, 72.34072430882337, 72.21263573835546, 72.22114559632422, 72.25515075157183, 72.15039144579399, 72.22911368388614, 72.3119562190233, 72.35206281480308, 72.31505249554667, 72.31474549154234, 72.27508746119776, 72.29176338828411, 72.16153243833875, 71.99881735047215, 72.37278040246152, 72.32137724551353, 72.28990024359747, 72.2099353247712, 72.24889401824284, 72.14848700307354, 72.1274604903228, 72.15406099641362, 72.16889619645852, 71.85926394185189, 71.26921398292495, 72.09686925059687, 71.41203713431179, 71.18770906971459, 71.26317269643839, 71.46594544995165, 71.66394619160673, 71.55258782314299, 71.68885820484513, 71.65245784749825, 71.64472999082103, 71.66792297963464, 71.44340715308148, 71.6870757194831, 71.67412462571818, 71.67111761057423, 71.6947500034728, 72.23199703269547, 72.32763764121421, 72.41673136694249, 72.30699543115739, 72.58628785989198, 72.37560790773551, 72.62775155836519, 72.61978463962099, 72.58558524653634, 72.54249873747199, 72.66195203675682], [479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 624.5753, 624.5753, 479.29004, 624.5753, 624.5753, 624.5753, 479.29004, 624.5753, 624.5753, 479.29004, 624.5753, 624.5753, 624.5753, 479.29004, 624.5753, 479.29004, 624.5753, 479.29004, 624.5753, 479.29004, 479.29004, 624.5753, 624.5753, 624.5753, 479.29004, 479.29004, 624.5753, 624.5753, 479.29004, 624.5753, 479.29004, 624.5753, 624.5753, 479.29004, 624.5753, 624.5753, 624.5753, 479.29004, 624.5753, 479.29004, 479.29004, 479.29004, 624.5753, 624.5753, 479.29004, 479.29004, 479.29004, 624.5753, 624.5753, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 724.99976, 724.99976, 724.99976, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 724.99976, 724.99976, 724.99976, 724.99976, 724.99976, 724.99976, 724.99976, 724.99976, 724.99976, 479.29004, 724.99976, 724.99976, 724.99976, 724.99976, 724.99976, 724.99976, 724.99976, 724.99976, 724.99976, 724.99976, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 724.99976, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 280.1672, 280.1672, 479.29004, 280.1672, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 280.1672, 479.29004, 479.29004, 280.1672, 280.1672, 280.1672, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 724.99976, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 479.29004, 724.99976, 724.99976, 479.29004, 724.99976, 724.99976, 724.99976, 360.89536, 360.89536, 360.89536, 360.89536, 724.99976, 724.99976, 360.89536, 360.89536, 724.99976, 360.89536, 360.89536, 360.89536, 479.29004, 479.29004, 280.1672, 479.29004, 280.1672, 280.1672, 280.1672, 280.1672, 280.1672, 280.1672, 280.1672], [3840.0, 3840.0, 3840.0, 3840.0, 3840.0, 3720.0, 3720.0, 3840.0, 3840.0, 3840.0, 3720.0, 3840.0, 3840.0, 3780.0, 3840.0, 3780.0, 3840.0, 3840.0, 3840.0, 3720.0, 3840.0, 3780.0, 3840.0, 3780.0, 3840.0, 3840.0, 3780.0, 3720.0, 3840.0, 3840.0, 3840.0, 3780.0, 3720.0, 3840.0, 3720.0, 3840.0, 3720.0, 3720.0, 3840.0, 3720.0, 3840.0, 3840.0, 3840.0, 3780.0, 3840.0, 3840.0, 3840.0, 3780.0, 3840.0, 3840.0, 3840.0, 3840.0, 3840.0, 3780.0, 3840.0, 3840.0, 3840.0, 3840.0, 3840.0, 3840.0, 3900.0, 3900.0, 3900.0, 3900.0, 3900.0, 3900.0, 3900.0, 3900.0, 3900.0, 3900.0, 3900.0, 3900.0, 3900.0, 3900.0, 3900.0, 3900.0, 3900.0, 3660.0, 3660.0, 3660.0, 3900.0, 3900.0, 3900.0, 3900.0, 3900.0, 3900.0, 3840.0, 3780.0, 3660.0, 3660.0, 3720.0, 3780.0, 3660.0, 3720.0, 3720.0, 3660.0, 3900.0, 3660.0, 3720.0, 3660.0, 3720.0, 3660.0, 3780.0, 3660.0, 3660.0, 3660.0, 3660.0, 3840.0, 3840.0, 3840.0, 3900.0, 3900.0, 3900.0, 3960.0, 3840.0, 4020.0, 3900.0, 3960.0, 4020.0, 3960.0, 4080.0, 4080.0, 3960.0, 4080.0, 3960.0, 4080.0, 3900.0, 4080.0, 4080.0, 4080.0, 3900.0, 3900.0, 4080.0, 4080.0, 4080.0, 3960.0, 3960.0, 3960.0, 4020.0, 3960.0, 4080.0, 3900.0, 3660.0, 4080.0, 4080.0, 3960.0, 3960.0, 4080.0, 3960.0, 3900.0, 4080.0, 4080.0, 3900.0, 3780.0, 3900.0, 3780.0, 3780.0, 3720.0, 4020.0, 4020.0, 4020.0, 4020.0, 3900.0, 3840.0, 4020.0, 4020.0, 3900.0, 4020.0, 4020.0, 4020.0, 4020.0, 4080.0, 4080.0, 3960.0, 4080.0, 4080.0, 4080.0, 4080.0, 4080.0, 4080.0, 4080.0], [87.63560920082658, 87.63560920082658, 87.63560920082658, 87.63560920082658, 87.63560920082658, 86.2554346113913, 86.2554346113913, 87.63560920082658, 87.63560920082658, 87.63560920082658, 86.2554346113913, 87.63560920082658, 87.63560920082658, 86.94826047713663, 87.63560920082658, 86.94826047713663, 87.63560920082658, 87.63560920082658, 87.63560920082658, 86.2554346113913, 87.63560920082658, 86.94826047713663, 87.63560920082658, 86.94826047713663, 87.63560920082658, 87.63560920082658, 86.94826047713663, 86.2554346113913, 87.63560920082658, 87.63560920082658, 87.63560920082658, 86.94826047713663, 86.2554346113913, 87.63560920082658, 86.2554346113913, 87.63560920082658, 86.2554346113913, 86.2554346113913, 87.63560920082658, 86.2554346113913, 87.63560920082658, 87.63560920082658, 87.63560920082658, 86.94826047713663, 87.63560920082658, 87.63560920082658, 87.63560920082658, 86.94826047713663, 87.63560920082658, 87.63560920082658, 87.63560920082658, 87.63560920082658, 87.63560920082658, 86.94826047713663, 87.63560920082658, 87.63560920082658, 87.63560920082658, 87.63560920082658, 87.63560920082658, 87.63560920082658, 88.31760866327846, 88.31760866327846, 88.31760866327846, 88.31760866327846, 88.31760866327846, 88.31760866327846, 88.31760866327846, 88.31760866327846, 88.31760866327846, 88.31760866327846, 88.31760866327846, 88.31760866327846, 88.31760866327846, 88.31760866327846, 88.31760866327846, 88.31760866327846, 88.31760866327846, 85.55699854482975, 85.55699854482975, 85.55699854482975, 88.31760866327846, 88.31760866327846, 88.31760866327846, 88.31760866327846, 88.31760866327846, 88.31760866327846, 87.63560920082658, 86.94826047713663, 85.55699854482975, 85.55699854482975, 86.2554346113913, 86.94826047713663, 85.55699854482975, 86.2554346113913, 86.2554346113913, 85.55699854482975, 88.31760866327846, 85.55699854482975, 86.2554346113913, 85.55699854482975, 86.2554346113913, 85.55699854482975, 86.94826047713663, 85.55699854482975, 85.55699854482975, 85.55699854482975, 85.55699854482975, 87.63560920082658, 87.63560920082658, 87.63560920082658, 88.31760866327846, 88.31760866327846, 88.31760866327846, 88.99438184514796, 87.63560920082658, 89.66604708583958, 88.31760866327846, 88.99438184514796, 89.66604708583958, 88.99438184514796, 90.33271832508971, 90.33271832508971, 88.99438184514796, 90.33271832508971, 88.99438184514796, 90.33271832508971, 88.31760866327846, 90.33271832508971, 90.33271832508971, 90.33271832508971, 88.31760866327846, 88.31760866327846, 90.33271832508971, 90.33271832508971, 90.33271832508971, 88.99438184514796, 88.99438184514796, 88.99438184514796, 89.66604708583958, 88.99438184514796, 90.33271832508971, 88.31760866327846, 85.55699854482975, 90.33271832508971, 90.33271832508971, 88.99438184514796, 88.99438184514796, 90.33271832508971, 88.99438184514796, 88.31760866327846, 90.33271832508971, 90.33271832508971, 88.31760866327846, 86.94826047713663, 88.31760866327846, 86.94826047713663, 86.94826047713663, 86.2554346113913, 89.66604708583958, 89.66604708583958, 89.66604708583958, 89.66604708583958, 88.31760866327846, 87.63560920082658, 89.66604708583958, 89.66604708583958, 88.31760866327846, 89.66604708583958, 89.66604708583958, 89.66604708583958, 89.66604708583958, 90.33271832508971, 90.33271832508971, 88.99438184514796, 90.33271832508971, 90.33271832508971, 90.33271832508971, 90.33271832508971, 90.33271832508971, 90.33271832508971, 90.33271832508971]]

gridx = x_y_z_t_sigma[0]
gridy = x_y_z_t_sigma[1]
sumt = x_y_z_t_sigma[3]

mx = [i*sstep for i in gridx]
my = [i*sstep for i in gridy]
mz = x_y_z_t_sigma[2] #z initunit is m

Rt_list = [] #drop time 
pi = 3.14
for i in range(1,len(mx)):
    a = max(mx[:i+1])-min(mx[:i+1])
    b = max(my[:i+1])-min(my[:i+1])
    Rt = math.sqrt(a*b/3.14)
    Rt_list.append(Rt)

# 排放时长te
te = [i+1 for i in range(len(Rt_list))] # 单位s

tenew = [i/10 for i in range(te[0]*10,te[-1]*10+1)]
# 下落时长td
Rdnew = [339.32 for i in range(len(tenew))] 
f=interpolate.interp1d(te,Rt_list,kind='slinear') # Rt 与te的插值函数
Rtnew = f(tenew)
print(len(tenew))
print(len(Rtnew))
print(len(Rdnew))
DDC1_list = []  # te = 60
DDC2_list = []  # te = 1
DDC3_list = []  # area add
for i in range(len(tenew)):
    DDC1 = DDC_func(Rdnew[i],Rtnew[i],tenew[i])
    DDC2 = DDC_func(Rdnew[i],Rtnew[i])
    DDC3 = DDC3_fun(Rdnew[i],Rtnew[i])
    DDC1_list.append(DDC1)
    DDC2_list.append(DDC2)
    DDC3_list.append(DDC3)

    
# #plot
fig=plt.figure(figsize=(8,6),dpi=300)#添加画布
plt.text(1, 0.96, 'emit hight = 3000m,sedimentation velocity = 1.1m/s, drop time = 63.5min')
x = tenew
y1 = DDC1_list
y2 = DDC2_list
y3 = DDC3_list
plt.plot(x,y1,color = 'g',linewidth = 1.2,label='DDC1')
plt.plot(x,y2,color = 'r',linewidth = 1.2,label='DDC2')
plt.plot(x,y3,color = 'b',linewidth = 1.2,label='DDC3')

plt.xlim(0, 60)
plt.xticks([i*10 for i in range(7)],[i*10 for i in range(7)])
plt.ylim(0,1)
plt.yticks([0,0.25,0.50,0.75,1],[0,25,50,75,100])
plt.xlabel('te(min)',size = 13)
plt.ylabel('ratio(%)', size = 13)
plt.legend(loc="lower right")

plt.savefig(FigurePath + FigureName)



