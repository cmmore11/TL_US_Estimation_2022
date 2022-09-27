# -*- coding: utf-8 -*-
"""
This code was developed in order to take a given .npy or .csv file of Npts width and 
length of Mmean values obtained using mean_mag_finder.py and main_code_web.py previously
and using them to calculated the TL at each point and plot it.

You must have the following modules installed: pyshp

This code yields a plot of the clipped lat/long points that each TL value was calculated at,
 a contour plot of the clipped TL values, and a surface plot of the clipped TL values.
@author: cmmore11
"""
##############################################################################
##               Importing Modules, reading Given .py File                  ##
##############################################################################
import numpy as np    
#import matplotlib                             #general calculations
import matplotlib.pyplot as plt                    #general plotting
import shapefile                                   #needs pyshp package, used to read USA shapefile
from shapely.geometry.polygon import Polygon 
#from matplotlib.patches import Polygon as Polygonn             #used to clip points with shapefile 
from shapely.geometry import Point  
import matplotlib.patches as mpatches
import matplotlib.path as mpath
#from scipy.interpolate import interp1d
from matplotlib.colors import TwoSlopeNorm
plt.rcParams["font.family"] = "Times New Roman"
#matplotlib.use('Qt5Agg')
Mmean4 = np.genfromtxt('modal_2s_HI.csv',delimiter=",")
Mmean4[0,] = 1
MmeanW = np.genfromtxt ('modal_1s_HI.csv', delimiter=",")
MmeanW[0,] = 1
##############################################################################
##                        calculate TL at each point                        ##
##############################################################################
# make a grid of latitude-longitude values, these are based off of the web scraping coords
# make a grid of latitude-longitude values, these are based off of the web scraping coords
xmin, xmax, ymin, ymax = -160.3, -154.8, 18.9, 22.3 #DO NOT CHANGE
Nptsx = 110 #grid length/height
Nptsy = 68 
xx, yy = np.meshgrid(np.linspace(xmin,xmax,Nptsx), np.linspace(ymin,ymax,Nptsy))
positions = np.vstack([xx.ravel(), yy.ravel()])
xc = xx.flatten() #x-coord vector
yc = yy.flatten() #y-coord vector

#initialization of matrices
M_mean = np.zeros([Nptsx*Nptsy,1])
M_mean4 = np.zeros([Nptsx*Nptsy,1])
TL = np.zeros([Nptsx*Nptsy,])
TL4 = np.zeros([Nptsx*Nptsy,])
TLC = np.zeros([Nptsx*Nptsy,])
sdrop_HI = 20 #bars  
beta_HI = 3.8 #km/s
#Boundary of all Hawaiian Islands, Starting with east-most island
h1 = [(-155.897,20.241),(-155.842,20.2698),(-155.7678,20.2453),(-155.733,20.220),(-155.73,20.2048),(-155.649,20.162),(-155.633,20.1426),(-155.6107,20.1345),(-155.589,20.118),(-155.556,20.1287),(-155.4386,20.0908),(-155.2697,20.015357),(-155.129,19.906),(-155.102,19.8669),(-155.087717,19.8598),(-155.08,19.84688),(-155.0898,19.80748),(-155.0918,19.7370),(-155.085,19.726),(-155.0689,19.7244),(-155.0383,19.739588),(-155.0213,19.733205),(-155.005,19.7379),(-154.9766,19.6854),(-154.985,19.6598),(-154.98,19.63267),(-154.9434,19.6239),(-154.9469,19.607),(-154.926785,19.5927),(-154.8686,19.5507),(-154.8255,19.53477),(-154.8066,19.52086),(-154.824,19.4735),(-154.968,19.3557),(-155.05,19.324),(-155.0748,19.31095),(-155.133,19.273356),(-155.1998,19.257799),(-155.2664,19.272),(-155.337,19.2338),(-155.35,19.209),(-155.42296,19.1806),(-155.45,19.147),(-155.5,19.1339),(-155.556,19.0807),(-155.556,19.052),(-155.566,19.026),(-155.58,19.021),(-155.596,18.9704),(-155.619,18.9665),(-155.6811,18.9106),(-155.685,18.935),(-155.758,18.9768),(-155.8816,19.0314),(-155.921566,19.115777),(-155.88575,19.309),(-155.91185,19.3945),(-155.9105,19.4152),(-155.9215,19.4327),(-155.9235,19.4806),(-155.948,19.4858),(-155.99,19.6301),(-156.028,19.6498),(-156.029,19.6706),(-156.062,19.7307),(-156.049,19.787577),(-155.9787,19.84636),(-155.953999,19.85282),(-155.936,19.850237),(-155.855,19.965),(-155.829,19.97822),(-155.824,20.0269),(-155.8819,20.103045),(-155.899,20.154622),(-155.9,20.218)]
h2 = [(-156.4768,20.8944),(-156.38,20.92),(-156.329,20.945),(-156.284,20.946),(-156.27,20.9348),(-156.2537,20.937),(-156.2296,20.9267),(-156.224,20.91388),(-156.162,20.863),(-156.136,20.8638),(-156.11287,20.824698),(-156.087,20.8272),(-156.059,20.80994),(-155.989,20.781050),(-155.98,20.7628),(-155.9867,20.722138),(-156.008,20.6812),(-156.024,20.6822),(-156.042,20.6662),(-156.0491,20.6512),(-156.111,20.6403),(-156.135,20.621),(-156.147,20.6267),(-156.171,20.621998),(-156.191,20.6273),(-156.296,20.584),(-156.3356,20.5852),(-156.371,20.5752),(-156.410,20.58199),(-156.42,20.598),(-156.429,20.5955),(-156.44,20.6067),(-156.44,20.618),(-156.4546,20.6362),(-156.44,20.6514),(-156.4488,20.7208),(-156.4628,20.7778),(-156.488,20.796),(-156.5,20.7944),(-156.522,20.77723),(-156.54,20.7759),(-156.5858,20.79536),(-156.6105,20.81093),(-156.6236,20.80868),(-156.68337,20.8759),(-156.687,20.903),(-156.696,20.91599),(-156.69,20.9489),(-156.6667,21.005),(-156.64,21.01367),(-156.639,21.0238),(-156.611,21.02328),(-156.599,21.0311),(-156.579,21.01487),(-156.55954,21.0115),(-156.537,20.9852),(-156.529,20.9855)]
h3 = [(-156.554733,20.534494),(-156.531215,20.534173),(-156.5456,20.51279),(-156.555,20.51987),(-156.58598,20.50813),(-156.598,20.5252),(-156.607,20.516),(-156.643,20.5054),(-156.547,20.51697),(-156.6586,20.50202),(-156.6737,20.5023),(-156.6998,20.526998),(-156.674,20.556195),(-156.647,20.5618),(-156.587,20.601996),(-156.5685,20.6047),(-156.5525,20.595087),(-156.54065,20.5705),(-156.5468,20.567446)]
h4 = [(-157.058,20.9108),(-157.025,20.9269),(-156.9365,20.9243),(-156.8943,20.9127),(-156.8181,20.8435),(-156.805,20.821),(-156.807,20.80625),(-156.84,20.757779),(-156.889,20.7385),(-156.964,20.7324),(-156.9818,20.7552),(-156.9906,20.78635),(-156.985,20.8182),(-157.008,20.853057),(-157.05,20.87),(-157.0599,20.886)]
h5 = [(-157.31098,21.098),(-157.2897,21.144),(-157.251,21.1823),(-157.249,21.205997),(-157.26154,21.2169),(-157.2554,21.2239),(-157.2,21.2175),(-157.18775,21.2041),(-156.995,21.181),(-156.985,21.1874),(-156.973,21.2134),(-156.9629,21.2124),(-156.944,21.1743),(-156.9248,21.1663),(-156.915,21.1714),(-156.903,21.1615),(-156.886,21.1685),(-156.876,21.1602),(-156.804,21.17623),(-156.79,21.16998),(-156.7729,21.17735),(-156.738,21.174146),(-156.7324,21.157177),(-156.709,21.15974),(-156.722,21.1325),(-156.7678,21.0902),(-156.804,21.0665),(-156.84845,21.05404),(-156.87797,21.048596),(-157.0517,21.0973),(-157.09049,21.101455),(-157.2498,21.08736)]
h6 = [(-158.278,21.5757),(-158.1238,21.58655),(-158.104,21.597),(-158.064,21.6402),(-158.065,21.6475),(-158.0167,21.698563),(-157.981375,21.7113),(-157.9628,21.70622),(-157.916,21.6459),(-157.92095,21.6293),(-157.8864,21.5813),(-157.875,21.5725),(-157.8745,21.55598),(-157.864,21.5613),(-157.849,21.55734),(-157.8344,21.5265),(-157.8366,21.50943),(-157.852,21.504),(-157.8269,21.4567),(-157.8149,21.4532),(-157.7833,21.41326),(-157.7696,21.4158),(-157.763,21.4379),(-157.78,21.4471),(-157.765,21.4602),(-157.744,21.4567),(-157.7226,21.4606),(-157.7442,21.419),(-157.737,21.4037),(-157.7156,21.39416),(-157.7064,21.38036),(-157.707,21.3557),(-157.69199,21.331117),(-157.65,21.3084),(-157.657,21.29338),(-157.68667,21.2758),(-157.692,21.26458),(-157.711,21.2607),(-157.715,21.2819),(-157.789,21.25866),(-157.79,21.255),(-157.81999,21.2552),(-157.825,21.27458),(-157.835,21.277),(-157.874,21.299453),(-157.907427,21.29768),(-157.95858,21.314953),(-157.97,21.3233),(-157.97777,21.3158),(-158.108,21.294),(-158.1366,21.3696),(-158.17792,21.4027),(-158.189,21.444581),(-158.2172,21.46583),(-158.22015,21.47637),(-158.2318,21.4832),(-158.231133,21.5366)]
h7 = [(-159.7798,22.0662),(-159.7468,22.0955),(-159.7276,22.15022),(-159.63697,22.1871),(-159.583,22.2278),(-159.544,22.22275),(-159.5412,22.217),(-159.523,22.217),(-159.508,22.2034),(-159.498,22.2105),(-159.498,22.2223),(-159.484,22.232),(-159.436,22.22377),(-159.4245,22.22329),(-159.402,22.2306),(-159.3858,22.2198),(-159.3501,22.21296),(-159.313387,22.1845),(-159.3046,22.15765),(-159.3057,22.1484),(-159.2938,22.14684),(-159.296655,22.10557),(-159.321,22.0553),(-159.3339,22.0483),(-159.338787,22.0234),(-159.3319,21.95958),(-159.35063,21.955759),(-159.3484,21.937289),(-159.4459,21.86896),(-159.5043,21.8879),(-159.5334,21.885),(-159.57566,21.893326),(-159.5925,21.9024),(-159.60347,21.892529),(-159.67,21.953),(-159.762,21.9801),(-159.7747,22.006),(-159.78887,22.0161)]
h8 = [(-160.2034,21.779),(-160.2302,21.79224),(-160.246,21.8158),(-160.2459,21.8467),(-160.233,21.8617),(-160.23,21.8853),(-160.1859,21.925750),(-160.16698,21.931483),(-160.1613,21.938967),(-160.11874,21.96253),(-160.11,21.99485),(-160.09145,22.00487),(-160.06667,21.99898),(-160.05,21.9864),(-160.054,21.97336),(-160.065,21.973357),(-160.0826,21.92909),(-160.073369,21.8931),(-160.1529,21.863918),(-160.19414,21.80693),(-160.19243,21.80422)]

h1_big = Polygon(h1)
h2_maui = Polygon(h2)
h3_kahoolawe = Polygon(h3)
h4_lanai = Polygon(h4)
h5_molokai = Polygon(h5)
h6_oahu = Polygon(h6)
h7_kauai = Polygon(h7)
h8_nihau = Polygon(h8)

#main loop to calculate TL at each point
for i in range(6450,Nptsx*Nptsy-400):
    if MmeanW[i,] ==1:
        M_mean[i,], M_mean4[i,] = 7, 7
    else:
        M_mean[i,] = MmeanW[i,]
        M_mean4[i,] = Mmean4[i,]
    TL[i,] = min(1/((4.9*(10**6))*beta_HI*((sdrop_HI/(10**(1.5*(10.7+M_mean[i,]))))**(1/3))),16)
    TL4[i,] = min(1/((4.9*(10**6))*beta_HI*((sdrop_HI/(10**(1.5*(10.7+M_mean4[i,]))))**(1/3))),16)
    if h1_big.contains(Point(xc[i,],yc[i,])) == True:
         TLC[i,] = 12;
    elif h2_maui.contains(Point(xc[i,],yc[i,])) == True:
         TLC[i,] = 6;
    elif h3_kahoolawe.contains(Point(xc[i,],yc[i,])) == True:
         TLC[i,] = 6;
    elif h4_lanai.contains(Point(xc[i,],yc[i,])) == True:
         TLC[i,] = 6;
    elif h5_molokai.contains(Point(xc[i,],yc[i,])) == True:
         TLC[i,] = 6;
    elif h6_oahu.contains(Point(xc[i,],yc[i,])) == True:
         TLC[i,] = 4;
    elif h7_kauai.contains(Point(xc[i,],yc[i,])) == True:
         TLC[i,] = 4;   
    elif h8_nihau.contains(Point(xc[i,],yc[i,])) == True:
         TLC[i,] = 4;
    else:
         TLC[i,] = 8

    print('xi = %2f'%xc[i,])
    print('yi = %2f'%yc[i,])
    print('TLi = %2f'%(TL[i,]))
    print('iteration %2d complete'%(i+1))

##############################################################################
##                           Beginning Figure                               ##
##############################################################################
plt.figure(num=1,clear=True,figsize=[22.4,16])
ax = plt.axes() # add the axes

##############################################################################
##                             Making patches                               ##
##############################################################################
#making patch in shape of US by making a "polygon" with a lot of points
p1, p2, p3, p4 = plt.Polygon(h1), plt.Polygon(h2), plt.Polygon(h3), plt.Polygon(h4)
all_polys = [p1,p2,p3,p4]
#h1
poly_codes1 = [mpath.Path.MOVETO] + (len(h1) - 2) * [mpath.Path.LINETO] + [mpath.Path.CLOSEPOLY]
# create a Path from the polygon vertices
path1 = mpath.Path(h1, poly_codes1)
# create a Patch from the path
patch1 = mpatches.PathPatch(path1, facecolor='none', edgecolor='none')
h2
poly_codes2 = [mpath.Path.MOVETO] + (len(h2) - 2) * [mpath.Path.LINETO] + [mpath.Path.CLOSEPOLY]
path2 = mpath.Path(h2, poly_codes2)
patch2 = mpatches.PathPatch(path2, facecolor='none', edgecolor='none')
#h3
poly_codes3 = [mpath.Path.MOVETO] + (len(h3) - 2) * [mpath.Path.LINETO] + [mpath.Path.CLOSEPOLY]
path3 = mpath.Path(h3, poly_codes3)
patch3 = mpatches.PathPatch(path3, facecolor='none', edgecolor='none')
#h4
poly_codes4 = [mpath.Path.MOVETO] + (len(h4) - 2) * [mpath.Path.LINETO] + [mpath.Path.CLOSEPOLY]
path4 = mpath.Path(h4, poly_codes4)
patch4 = mpatches.PathPatch(path4, facecolor='none', edgecolor='none')
#h5
poly_codes5 = [mpath.Path.MOVETO] + (len(h5) - 2) * [mpath.Path.LINETO] + [mpath.Path.CLOSEPOLY]
path5 = mpath.Path(h5, poly_codes5)
patch5 = mpatches.PathPatch(path5, facecolor='none', edgecolor='none')
#h3
poly_codes6 = [mpath.Path.MOVETO] + (len(h6) - 2) * [mpath.Path.LINETO] + [mpath.Path.CLOSEPOLY]
path6 = mpath.Path(h6, poly_codes6)
patch6 = mpatches.PathPatch(path6, facecolor='none', edgecolor='none')
#h3
poly_codes7 = [mpath.Path.MOVETO] + (len(h7) - 2) * [mpath.Path.LINETO] + [mpath.Path.CLOSEPOLY]
path7 = mpath.Path(h7, poly_codes7)
patch7 = mpatches.PathPatch(path7, facecolor='none', edgecolor='none')
#h3
poly_codes8 = [mpath.Path.MOVETO] + (len(h8) - 2) * [mpath.Path.LINETO] + [mpath.Path.CLOSEPOLY]
path8 = mpath.Path(h8, poly_codes8)
patch8 = mpatches.PathPatch(path8, facecolor='none', edgecolor='none')

plt.rc('axes', labelsize=45)
ax.tick_params(axis='both',which='major',width=2,length=10)

ax.add_patch(patch1)
ax.add_patch(patch2)
ax.add_patch(patch3)
ax.add_patch(patch4)
ax.add_patch(patch5)
ax.add_patch(patch6)
ax.add_patch(patch7)
ax.add_patch(patch8)

##############################################################################
##                          Making Arrays to Plot                           ##
##############################################################################
#need grid of values for contour plot AND surface plot later
xx2 = xc.reshape(Nptsy,Nptsx)
yy2 = yc.reshape(Nptsy,Nptsx)
TL2 = TL.reshape(Nptsy,Nptsx)
TL4 = TL4.reshape(Nptsy,Nptsx)
TLC2 = TLC.reshape(Nptsy,Nptsx)
M_mean = M_mean.reshape(Nptsy,Nptsx)
M_mean4 = M_mean4.reshape(Nptsy,Nptsx)

###############################################################################
################          Plot what you want         ##########################
#####comment/uncomment groups of lines to toggle on/off plotting that figure.##
#####make sure to change appropriate/corresponding figure title below #########
###############################################################################
# plt.figure(1)
# ax = plt.gca()

##Plotting modal magnitude difference at different spectral periods for Hawaii (Figure 10a)#################################################################
# tick_marks = [-6.7,-6,-4,-2,-1,-0.5,-0.2,0.2,0.5,1,1.25]
# norm_num = TwoSlopeNorm(vmin = -6.7, vcenter = 0, vmax = 6.7)
# cont =   ax.contourf(xx2,yy2,M_mean4-M_mean,tick_marks,norm=norm_num,cmap='coolwarm_r')
# cont2 = ax.contourf(xx2,yy2,M_mean4-M_mean,tick_marks,norm=norm_num,cmap='coolwarm_r')
# cont3 = ax.contourf(xx2,yy2,M_mean4-M_mean,tick_marks,norm=norm_num,cmap='coolwarm_r')
# cont4 = ax.contourf(xx2,yy2,M_mean4-M_mean,tick_marks,norm=norm_num,cmap='coolwarm_r')
# cont5 = ax.contourf(xx2,yy2,M_mean4-M_mean,tick_marks,norm=norm_num,cmap='coolwarm_r')
# cont6 = ax.contourf(xx2,yy2,M_mean4-M_mean,tick_marks,norm=norm_num,cmap='coolwarm_r')
# cont7 = ax.contourf(xx2,yy2,M_mean4-M_mean,tick_marks,norm=norm_num,cmap='coolwarm_r')
# cont8 = ax.contourf(xx2,yy2,M_mean4-M_mean,tick_marks,norm=norm_num,cmap='coolwarm_r')

##Plotting TL difference at different spectral periods for Hawaii (Figure 10b)#######################################################################
# tick_marks = [-16,-14,-12,-10,-8,-6,-4,-2,-0.5,0.5,2,4,6,8,10,11.05]
# norm_num = TwoSlopeNorm(vmin = -16, vcenter = 0, vmax = 16)
# cont =   ax.contourf(xx2, yy2, TL4-TL2,tick_marks,norm=norm_num, cmap='coolwarm_r')
# cont2 = ax.contourf(xx2, yy2, TL4-TL2,tick_marks,norm=norm_num, cmap='coolwarm_r')
# cont3 = ax.contourf(xx2, yy2, TL4-TL2,tick_marks,norm=norm_num, cmap='coolwarm_r')
# cont4 = ax.contourf(xx2, yy2, TL4-TL2,tick_marks,norm=norm_num, cmap='coolwarm_r')
# cont5 = ax.contourf(xx2, yy2, TL4-TL2,tick_marks,norm=norm_num, cmap='coolwarm_r')
# cont6 = ax.contourf(xx2, yy2, TL4-TL2,tick_marks,norm=norm_num, cmap='coolwarm_r')
# cont7 = ax.contourf(xx2, yy2, TL4-TL2,tick_marks,norm=norm_num, cmap='coolwarm_r')
# cont8 = ax.contourf(xx2, yy2, TL4-TL2,tick_marks,norm=norm_num, cmap='coolwarm_r')

##Plotting modal magnitude at spectral period of 1s for Hawaii (Figure 11)###########################################################################
# tick_marks = [6.1,6.25,6.5,6.75,7,7.25,7.5,7.7]
# cont  =  ax.contourf(xx2,yy2,M_mean,tick_marks,cmap='Purples')
# cont2 = ax.contourf(xx2,yy2,M_mean,tick_marks,cmap='Purples')
# cont3 = ax.contourf(xx2,yy2,M_mean,tick_marks,cmap='Purples')
# cont4 = ax.contourf(xx2,yy2,M_mean,tick_marks,cmap='Purples')
# cont5 = ax.contourf(xx2,yy2,M_mean,tick_marks,cmap='Purples')
# cont6 = ax.contourf(xx2,yy2,M_mean,tick_marks,cmap='Purples')
# cont7 = ax.contourf(xx2,yy2,M_mean,tick_marks,cmap='Purples')
# cont8 = ax.contourf(xx2,yy2,M_mean,tick_marks,cmap='Purples')


##Plotting final TL from this study for Hawaii (Figure 15)###########################################################################
tick_marks = [4,6,8,10,12,14,16]
norm_num = TwoSlopeNorm(vmin = 0, vcenter = 10, vmax = 20)
cont   = ax.contourf(xx2, yy2, TL2,tick_marks,norm =norm_num,cmap='Blues')
cont2 = ax.contourf(xx2, yy2, TL2,tick_marks,norm =norm_num,cmap='Blues')
cont3 = ax.contourf(xx2, yy2, TL2,tick_marks,norm =norm_num,cmap='Blues')
cont4 = ax.contourf(xx2, yy2, TL2,tick_marks,norm =norm_num,cmap='Blues')
cont5 = ax.contourf(xx2, yy2, TL2,tick_marks,norm =norm_num,cmap='Blues')
cont6 = ax.contourf(xx2, yy2, TL2,tick_marks,norm =norm_num,cmap='Blues')
cont7 = ax.contourf(xx2, yy2, TL2,tick_marks,norm =norm_num,cmap='Blues')
cont8 = ax.contourf(xx2, yy2, TL2,tick_marks,norm =norm_num,cmap='Blues')


##Plotting the difference between approx code and our TL (Figure 17)#############################################
# tick_marks = [-10.0,-9,-8,-7,-5,-3,-2,-1,-0.5,-0.25]
# norm_num = TwoSlopeNorm(vmin = -12, vcenter = 0, vmax = 12)
# cont   = ax.contourf(xx2, yy2, TLC2-TL2,tick_marks,norm=norm_num,cmap='coolwarm_r')
# cont2 = ax.contourf(xx2, yy2, TLC2-TL2,tick_marks,norm=norm_num,cmap='coolwarm_r')
# cont3 = ax.contourf(xx2, yy2, TLC2-TL2,tick_marks,norm=norm_num,cmap='coolwarm_r')
# cont4 = ax.contourf(xx2, yy2, TLC2-TL2,tick_marks,norm=norm_num,cmap='coolwarm_r')
# cont5 = ax.contourf(xx2, yy2, TLC2-TL2,tick_marks,norm=norm_num,cmap='coolwarm_r')
# cont6 = ax.contourf(xx2, yy2, TLC2-TL2,tick_marks,norm=norm_num,cmap='coolwarm_r')
# cont7 = ax.contourf(xx2, yy2, TLC2-TL2,tick_marks,norm=norm_num,cmap='coolwarm_r')
# cont8 = ax.contourf(xx2, yy2, TLC2-TL2,tick_marks,norm=norm_num,cmap='coolwarm_r')

# setting clip path##############################################################################################

for c in cont.collections:
    c.set_clip_path(patch1)
for c2 in cont2.collections:
    c2.set_clip_path(patch2)
for c3 in cont3.collections:
    c3.set_clip_path(patch3)
for c4 in cont4.collections:
    c4.set_clip_path(patch4)
for c5 in cont5.collections:
    c5.set_clip_path(patch5)
for c6 in cont6.collections:
    c6.set_clip_path(patch6)
for c7 in cont7.collections:
    c7.set_clip_path(patch7)
for c8 in cont8.collections:
    c8.set_clip_path(patch8)
    
plt.rc('axes', labelsize=45)
ax.tick_params(axis='both',which='major',width=2,length=10)

cbar = plt.colorbar(cont,fraction=0.02, pad=0.02)
cbar.ax.tick_params(axis='both',which='major',width=2,length=10)
#################################################################################
################uncomment 1 of the cbar.set_label lines########################
###############################################################################

##Plotting modal magnitude difference at different spectral periods for Hawaii (Figure 10a)#################################################################
#cbar.set_label("$M_w$$_($$_4$$_s$$_)$-$M_w$$_($$_1$$_s$$_)$",rotation=270,labelpad=50)

##Plotting TL difference at different spectral periods for Hawaii (Figure 10b)#######################################################################
#cbar.set_label("$T_L$$_($$_2$$_s$$_)$- $T_L$$_($$_1$$_s$$_)$", rotation=270,labelpad=50)

##Plotting modal magnitude at spectral period of 1s for Hawaii (Figure 11)###########################################################################
#cbar.set_label('Modal Magnitude, $M_w$',rotation=270,labelpad=60)

##Plotting final TL from this study for Hawaii (Figure 15)###########################################################################
cbar.set_label("$T_L$", rotation=270,labelpad=50)

##Plotting the difference between approx code and our TL (Figure 17)#############################################
#cbar.set_label("$T_L$$_($$_N$$_E$$_H$$_R$$_P$$_)$$_s$$_i$$_m$$_.$ - $T_L$", rotation=270,labelpad=45)

##############################################################################
##                          Plotting Shapefiles                           ##
##############################################################################
#Border of USA
sf = shapefile.Reader('cb_2018_us_nation_5m')
for shape in list(sf.iterShapes()):
    npoints=len(shape.points) # total points
    nparts = len(shape.parts) # total parts
    
    if nparts == 1:
        x_lon = np.zeros((len(shape.points),1))
        y_lat = np.zeros((len(shape.points),1))
        coord = np.zeros((len(shape.points),1))
        for ip in range(len(shape.points)):
            x_lon[ip] = shape.points[ip][0]
            y_lat[ip] = shape.points[ip][1] 
        plt.plot(x_lon,y_lat,color='k',linewidth=2)
    else: # loop over parts of each shape, plot separately
        for ip in range(nparts): # loop over parts, plot separately
            i0=shape.parts[ip]
            if ip < nparts-1:
              i1 = shape.parts[ip+1]-1
            else:
              i1 = npoints
            seg=shape.points[i0:i1+1]
            x_lon = np.zeros((len(seg),1))
            y_lat = np.zeros((len(seg),1))
            coord = np.zeros((len(seg),1))
            for ip in range(len(seg)):
              x_lon[ip] = seg[ip][0]
              y_lat[ip] = seg[ip][1] 
            plt.plot(x_lon,y_lat,color='k',linewidth=2.5)

plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.rc('axes', labelsize=40)
plt.rc('font',size=40)
csfont = {'fontname':'Times New Roman'}
ax.set_xlabel('Longitude',fontsize=45)
ax.set_ylabel('Latitude',fontsize=45)
plt.show()

#plot clipped points #optional
#t_l_points = ax.scatter(xc,yc,s=0.5,color='k')
#t_l_points.set_clip_path(patch)


