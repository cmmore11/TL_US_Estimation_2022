# -*- coding: utf-8 -*-
"""
This code was developed in order to take a given .csv file of length = (Nptsx X nptsy) of Modal magnitude values 
(MmeanW) and using it to calculate the TL at each point and plot it in the conterminous US

You must have the following modules installed: pyshp

This code yields a plot of the clipped lat/long points that each TL value was calculated at,
 and a contour plot of the clipped TL values.
@author: Christie Assadollahi, cmmore11@memphis.edu
"""
##############################################################################
##               Importing Modules, reading Given .py File                  ##
##############################################################################
import numpy as np    
#import matplotlib                             #general calculations
import matplotlib.pyplot as plt                    #general plotting
import shapefile                                   #needs pyshp package, used to read USA shapefile
from shapely.geometry.polygon import Polygon 
from matplotlib.patches import Polygon as Polygonn             #used to clip points with shapefile 
from shapely.geometry import Point  
import matplotlib.patches as mpatches
import matplotlib.path as mpath
from scipy.interpolate import interp1d
from matplotlib.colors import TwoSlopeNorm
plt.rcParams["font.family"] = "Times New Roman"
#matplotlib.use('Qt5Agg')
Mmean4 = np.genfromtxt('modal_4s_conterminous.csv',delimiter=",")
MmeanW = np.genfromtxt ('modal_2s_conterminous.csv', delimiter=",")
##############################################################################
##                        calculate TL at each point                        ##
##############################################################################
# make a grid of latitude-longitude values, these are based off of the web scraping coords
xmin, xmax, ymin, ymax = -125, -65, 24.6, 50 #DO NOT CHANGE
Nptsx = 120 #grid length/height
Nptsy = 52 
xx, yy = np.meshgrid(np.linspace(xmin,xmax,Nptsx), np.linspace(ymin,ymax,Nptsy))
positions = np.vstack([xx.ravel(), yy.ravel()])
xc = xx.flatten() #x-coord vector
yc = yy.flatten() #y-coord vector
# Saving the array in a text file
np.savetxt("xc.txt", xc)
np.savetxt("yc.txt",yc)
# file = open("xc.txt", "w+")
# content = str(xc)
# file.write(content)
# file.close()
#initialization of matrices
M_mean = np.zeros([Nptsx*Nptsy,1])
M_mean4 = np.zeros([Nptsx*Nptsy,1])
StressD = np.zeros([Nptsx*Nptsy,1])
StressD4 = np.zeros([Nptsx*Nptsy,1])
beta = np.zeros([Nptsx*Nptsy,1])
TL = np.zeros([Nptsx*Nptsy,])
TL4 = np.zeros([Nptsx*Nptsy,])
TLC = np.zeros([Nptsx*Nptsy,])


#function for stress drop in WUS
M_west = np.array([3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,9.5])
Stress_D_west = np.array([80,80,105,150,235,235,203,180,147,120,95,95]).squeeze()
fw = interp1d(M_west,Stress_D_west, kind = 'linear',axis=-1)#equation for stress drop from Zandieh, Pezeshk Campbell 2018 graph pg. 8

#function for stress drop in CEUS
M_cena = np.array([4,4.5,5,5.5,6,6.5,7,7.5,8,9.5])
Stress_D_CENA = np.array([102,115,122,133,128,132,126,121,115,115]).squeeze() #median
fc = interp1d(M_cena,Stress_D_CENA, kind = 'linear')

#function for Peterson's boundary between CEUS
bound_x = [-113.7,-113.35,-112.2,-110.9,-110.45,-110,-110.1,-109.9,-109.85,-110.28,-110.72,-110.7,-111,-111.5,-111.6,-111.4,-110.9,-110.9,-111.5,-111.35,-111,-110.25,-109.7,-109.3,-109,-108.45,-107.6,-107.15,-106.9,-106.85,-107,-107.55,-107.75,-107.8,-107.5,-106.25,-105.9,-105.62,-105.28,-105.09,-105.12,-105.17,-105.22,-105.2,-105.55,-105.55,-105.35,-104.5,-102.1,-101,-101] 
bound_y = [50,49,47.55,46.45,45.75,45.5,44.9,44.55,44.15,43.17,42.18,39.55,39.4,38.5,37.95,37.45,37.05,36.8,36.2,35.7,35.25,34.8,34.55,34.55,34.7,35,35.5,35.9,36.25,36.65,37,37.3,37.55,37.85,38.4,38.55,38.3,38.08,37.68,37.29,36.81,36.59,35.95,35,34.35,33.85,33.25,31.7,29.6,29.2,24]
bound_line = interp1d(bound_y,bound_x,kind = 'linear',axis=-1)

#Boundary of entire conterminous US
usapoly = Polygon([(-67.8,45.7),(-67,44.94),(-67,44.7),(-68.6,44.15),(-70.2,43.67),(-70.81,42.8),(-70.4,43.48),(-70.86,42.9),(-70.5,41.62),(-72.2,40.98),(-72.92,40.71),(-74.04,40.55),(-74.22,40.47),(-73.99,40.38),(-74.06,39.81),(-74.96,38.93),(-75.57,39.48),(-75.10,38.75),(-75.06,38.45),(-75.93,37.12),(-76.01,37.30),(-75.87,37.93),(-76,38.305),(-76,37.92),(-76.16,38.23),(-76.38,38.78),(-75.97,39.565),(-76.496,39.16),(-76.22,37.27),(-75.35,35.5),(-79.76,32.805),(-81.41,31.03),(-80.1,26.63),(-80,25.8),(-80.3,25.3),(-80.79,25.2),(-81.15,25.15),(-82.37,26.71),(-82.72,29.05),(-84.40,30.04),(-85,29.6),(-86.5,30.3),(-89.58,30.19),(-89.23,29.05),(-90.8,29.05),(-92,29.75),(-92.25,29.6),(-93,29.7),(-93.89,29.7),(-96.53,28.35),(-96.90,28.26),(-97.4,26.99),(-97.16,25.97),(-97.41,25.86),(-99.09,26.42),(-99.44,27.05),(-99.54,27.6),(-99.68,27.66),(-100.28,28.28),(-100.68,29.1),(-101.43,29.79),(-102.69,29.73),(-103.18,28.99),(-104,29.33),(-104.65,29.84),(-104.95,30.64),(-106.13,31.42),(-106.52,31.83),(-108.2,31.8),(-108.22,31.37),(-111.1,31.35),(-114.8,32.51),(-114.82,32.72),(-117.12,32.54),(-117.25,33.16),(-118.13,33.76),(-118.250,33.708),(-118.4,33.79),(-118.547,34.05),(-118.80,34.014),(-119.23,34.16),(-119.30,34.28),(-119.6,34.43),(-120.46,34.45),(-120.65,34.54),(-120.66,35.13),(-120.90,35.25),(-120.879,35.41),(-121.87,36.30),(-121.95,36.61),(-121.88,36.61),(-121.79,36.82),(-121.95,36.99),(-122.08,36.97),(-122.41,37.22),(-122.51,37.54),(-122.51,37.82),(-122.91,38.03),(-123.01,38),(-122.96,38.22),(-123.75,38.96),(-123.86,39.81),(-124.35,40.26),(-124.41,40.44),(-124.03,41.37),(-124.55,42.83),(-124.38,43.31),(-124.24,43.5),(-124.06,44.81),(-123.97,45.81),(-124.214,47.097),(-124.64,47.92),(-124.71,48.4),(-122.9,48.15),(-122.67,48.41),(-123.19,48.59),(-122.77,49),(-95.15,49),(-95.16,49.37),(-94.82,49.33),(-94.7,49.01),(-94.68,48.78),(-93.29,48.64),(-92.47,48.42),(-92.35,48.23),(-92.28,48.36),(-91.56,48.08),(-90.86,48.25),(-90.77,48.1),(-89.53,48.02),(-90.22,47.78),(-90.87,47.57),(-92.06,46.84),(-92.1,46.79),(-91.99,46.7),(-90.975,46.95),(-90.79,46.9),(-90.95,46.62),(-90.74,46.69),(-90.48,46.6),(-89.07,46.89),(-88.43,47.37),(-87.73,47.43),(-88.23,47.21),(-88.47,46.76),(-88.22,46.9),(-87.68,46.83),(-87.39,46.51),(-87.02,46.53),(-86.91,46.45),(-86.17,46.66),(-85.03,46.76),(-85.05,46.5),(-84.13,46.52),(-84.06,46.17),(-83.88,45.99),(-83.8,46.1),(-83.62,46.09),(-83.48,46),(-83.54,45.92),(-84.38,45.95),(-84.65,46.06),(-84.72,45.85),(-85.39,46.09),(-85.64,45.97),(-85.995,45.96),(-86.65,45.57),(-86.7,45.66),(-86.55,45.89),(-86.97,45.68),(-87.03,45.82),(-87.75,45),(-88.05,44.57),(-87.95,44.53),(-86.95,45.42),(-86.83,45.4),(-87.54,44.33),(-87.51,44.20),(-87.9,43.37),(-87.84,42.46),(-87.52,41.71),(-87.05,41.67),(-86.63,41.90),(-86.26,42.46),(-86.35,44.75),(-86.19,42.96),(-86.52,43.63),(-86.5,44.06),(-86.29,44.31),(-86.25,44.67),(-86.06,44.9),(-85.63,45.14),(-85.66,44.83),(-85.52,44.74),(-85.33,45.29),(-85.03,45.68),(-84.75,45.77),(-83.44,45.32),(-83.29,45.07),(-83.45,44.94),(-83.29,44.83),(-83.3,44.4),(-83.89,43.87),(-83.794,43.62),(-82.94,44.07),(-82.68,43.86),(-82.43,43),(-82.66,42.65),(-83.45,41.78),(-82.92,41.57),(-81.72,41.5),(-80.72,41.9), (-79.23,42.52),(-78.87,42.85),(-79,43.28),(-76.68,43.3),(-76.26,43.56),(-76.28,44.17),(-74.92,45.06),(-71.933,44.94),(-71.43,45.21),(-70.655,45.58),(-70,46.7),(-69.27,47.43),(-68.33,47.35),(-67.8,47.08),(-67.78,45.91)]) #manually enter some coords for rough outline of USA

#function for NEHRPs approximation/simplification of TL versus magnitude
M_nehrp = np.array([3.5,6.5,6.51,7,7.01,7.5,7.51,8,8.01,9.5])
Tc = np.array([4,4,6,6,8,8,12,12,16,16])
f_nehrp = interp1d(M_nehrp,Tc,kind = 'linear')
longlat = np.array([-85,40])
#main loop to calculate TL at each point
for i in range(0,Nptsx*Nptsy): #
    if Mmean4[i,] ==1:
        M_mean[i,] = 7
        M_mean4[i,] = 7
    else:
        M_mean[i,] = MmeanW[i,]
        M_mean4[i,] = Mmean4[i,]
    #M_mean[i,] = mean_mag_finder(yc[i,],xc[i,],760,2,2475)
    if usapoly.contains(Point(xc[i,],yc[i,])) == True:
        if xc[i,] < bound_line(yc[i,]):
            StressD[i,] =fw(M_mean[i,])
            StressD4[i,] = fw(M_mean4[i,])
            beta[i,] = 3.5
        else:
            StressD[i,] = fc(M_mean[i,])
            StressD4[i,] = fc(M_mean4[i,])
            beta[i,] = 3.7
        TL[i,] = min(1/((4.9*(10**6))*beta[i,]*((StressD[i,]/(10**(1.5*(10.7+M_mean[i,]))))**(1/3))),16)
        TL4[i,] = min(1/((4.9*(10**6))*beta[i,]*((StressD4[i,]/(10**(1.5*(10.7+M_mean4[i,]))))**(1/3))),16)
        TLC[i,] = min(f_nehrp(M_mean[i,]),16)
    else:
        TL[i,] = 0
        TL4[i,] = 0
        TLC[i,] = 0
    
        
    # TLa[i,] =  
     # # min(10**(-1.25+0.3*M_mean[i,]),16)
    print('xi = %2f'%xc[i,])
    print('yi = %2f'%yc[i,])
    print('TLi-TLCi = %2f'%(TL[i,]-TLC[i,]))
    print('iteration %2d complete'%(i+1))

##############################################################################
##                           Beginning Figure                               ##
##############################################################################
plt.figure(num=1,clear=True,figsize=[22.4,16])
ax = plt.axes() # add the axes
# ax.set_aspect(1.4)

##############################################################################
##                Making and Clipping using USA Polygon                     ##
##############################################################################
#making patch in shape of US by making a "polygon" with a lot of points
#coords from Bing maps, start at Maine and go clockwise around boundary
coords = [(-67.8,45.7),(-67,44.94),(-67,44.7),(-68.6,44.15),(-70.2,43.67),(-70.81,42.8),(-70.4,43.48),(-70.86,42.9),(-70.5,41.62),(-72.2,40.98),(-72.92,40.71),(-74.04,40.55),(-74.22,40.47),(-73.99,40.38),(-74.06,39.81),(-74.96,38.93),(-75.57,39.48),(-75.10,38.75),(-75.06,38.45),(-75.93,37.12),(-76.01,37.30),(-75.87,37.93),(-76,38.305),(-76,37.92),(-76.16,38.23),(-76.38,38.78),(-75.97,39.565),(-76.496,39.16),(-76.22,37.27),(-75.35,35.5),(-79.76,32.805),(-81.41,31.03),(-80.1,26.63),(-80,25.8),(-80.3,25.3),(-80.79,25.2),(-81.15,25.15),(-82.37,26.71),(-82.72,29.05),(-84.40,30.04),(-85,29.6),(-86.5,30.3),(-89.58,30.19),(-89.23,29.05),(-90.8,29.05),(-92,29.75),(-92.25,29.6),(-93,29.7),(-93.89,29.7),(-96.53,28.35),(-96.90,28.26),(-97.4,26.99),(-97.16,25.97),(-97.41,25.86),(-99.09,26.42),(-99.44,27.05),(-99.54,27.6),(-99.68,27.66),(-100.28,28.28),(-100.68,29.1),(-101.43,29.79),(-102.69,29.73),(-103.18,28.99),(-104,29.33),(-104.65,29.84),(-104.95,30.64),(-106.13,31.42),(-106.52,31.83),(-108.2,31.8),(-108.22,31.37),(-111.1,31.35),(-114.8,32.51),(-114.82,32.72),(-117.12,32.54),(-117.25,33.16),(-118.13,33.76),(-118.250,33.708),(-118.4,33.79),(-118.547,34.05),(-118.80,34.014),(-119.23,34.16),(-119.30,34.28),(-119.6,34.43),(-120.46,34.45),(-120.65,34.54),(-120.66,35.13),(-120.90,35.25),(-120.879,35.41),(-121.87,36.30),(-121.95,36.61),(-121.88,36.61),(-121.79,36.82),(-121.95,36.99),(-122.08,36.97),(-122.41,37.22),(-122.51,37.54),(-122.51,37.82),(-122.91,38.03),(-123.01,38),(-122.96,38.22),(-123.75,38.96),(-123.86,39.81),(-124.35,40.26),(-124.41,40.44),(-124.03,41.37),(-124.55,42.83),(-124.38,43.31),(-124.24,43.5),(-124.06,44.81),(-123.97,45.81),(-124.214,47.097),(-124.64,47.92),(-124.71,48.4),(-122.9,48.15),(-122.67,48.41),(-123.19,48.59),(-122.77,49),(-95.15,49),(-95.16,49.37),(-94.82,49.33),(-94.7,49.01),(-94.68,48.78),(-93.29,48.64),(-92.47,48.42),(-92.35,48.23),(-92.28,48.36),(-91.56,48.08),(-90.86,48.25),(-90.77,48.1),(-89.53,48.02),(-90.22,47.78),(-90.87,47.57),(-92.06,46.84),(-92.1,46.79),(-91.99,46.7),(-90.975,46.95),(-90.79,46.9),(-90.95,46.62),(-90.74,46.69),(-90.48,46.6),(-89.07,46.89),(-88.43,47.37),(-87.73,47.43),(-88.23,47.21),(-88.47,46.76),(-88.22,46.9),(-87.68,46.83),(-87.39,46.51),(-87.02,46.53),(-86.91,46.45),(-86.17,46.66),(-85.03,46.76),(-85.05,46.5),(-84.13,46.52),(-84.06,46.17),(-83.88,45.99),(-83.8,46.1),(-83.62,46.09),(-83.48,46),(-83.54,45.92),(-84.38,45.95),(-84.65,46.06),(-84.72,45.85),(-85.39,46.09),(-85.64,45.97),(-85.995,45.96),(-86.65,45.57),(-86.7,45.66),(-86.55,45.89),(-86.97,45.68),(-87.03,45.82),(-87.75,45),(-88.05,44.57),(-87.95,44.53),(-86.95,45.42),(-86.83,45.4),(-87.54,44.33),(-87.51,44.20),(-87.9,43.37),(-87.84,42.46),(-87.52,41.71),(-87.05,41.67),(-86.63,41.90),(-86.26,42.46),(-86.35,44.75),(-86.19,42.96),(-86.52,43.63),(-86.5,44.06),(-86.29,44.31),(-86.25,44.67),(-86.06,44.9),(-85.63,45.14),(-85.66,44.83),(-85.52,44.74),(-85.33,45.29),(-85.03,45.68),(-84.75,45.77),(-83.44,45.32),(-83.29,45.07),(-83.45,44.94),(-83.29,44.83),(-83.3,44.4),(-83.89,43.87),(-83.794,43.62),(-82.94,44.07),(-82.68,43.86),(-82.43,43),(-82.66,42.65),(-83.45,41.78),(-82.92,41.57),(-81.72,41.5),(-80.72,41.9), (-79.23,42.52),(-78.87,42.85),(-79,43.28),(-76.68,43.3),(-76.26,43.56),(-76.28,44.17),(-74.92,45.06),(-71.933,44.94),(-71.43,45.21),(-70.655,45.58),(-70,46.7),(-69.27,47.43),(-68.33,47.35),(-67.8,47.08),(-67.78,45.91)] #manually enter some coords for rough outline of USA
polygons = Polygonn(coords)
poly_codes = [mpath.Path.MOVETO] + (len(coords) - 2) * [mpath.Path.LINETO] + [mpath.Path.CLOSEPOLY]
# create a Path from the polygon vertices
path = mpath.Path(coords, poly_codes)
# create a Patch from the path
patch = mpatches.PathPatch(path, facecolor='none', edgecolor='none')
plt.figure(1)
ax = plt.gca()
# add the patch to the axes
ax.add_patch(patch)

##############################################################################
##                          Making Arrays to Plot                           ##
##############################################################################
#need grid of values for contour plot AND surface plot later
xx2 = xc.reshape(Nptsy,Nptsx)
yy2 = yc.reshape(Nptsy,Nptsx)
TL2 = TL.reshape(Nptsy,Nptsx)
TL4 = TL4.reshape(Nptsy,Nptsx)
TLC2 = TLC.reshape(Nptsy,Nptsx)
StressD = StressD.reshape(Nptsy,Nptsx)
M_mean = M_mean.reshape(Nptsy,Nptsx)
M_mean4 = M_mean4.reshape(Nptsy,Nptsx)

###############################################################################
################          Plot what you want         ##########################
###############################################################################
# maxdif = np.max(M_mean)
# mindif = np.min(M_mean)
# print('max = %2f'%maxdif)
# print('min = %2f'%mindif)

##Plotting modal magnitude difference at different spectral periods for CONUS (Figure 8a)############################################################
#cont = ax.contourf(xx2,yy2,M_mean4-M_mean,[-0.25,-0.2,0.2,0.5,1,2,3,3.05],norm=TwoSlopeNorm(vmin = -3.05, vcenter = 0, vmax = 3.05),cmap='coolwarm_r')

##Plotting TL difference at different spectral periods for CONUS (Figure 8b)#########################################################################
#cont = ax.contourf(xx2, yy2, TL4-TL2,[-3,-2,-1,-0.5,0.5,1,2,3,5,7,9,11,13,13.5],norm=TwoSlopeNorm(vmin = -13.5, vcenter = 0, vmax = 13.5), cmap='coolwarm_r')

##Plotting modal magnitude at spectral period of 2s for CONUS (Figure 9)###########################################################################
#cont = ax.contourf(xx2,yy2,M_mean,[5.3,6,6.5,7,7.5,8,8.5,9,9.35],cmap='Purples')

##Plotting stress drop for CONUS (Figure 12) ########################################################################################################
#cont = ax.contourf(xx2,yy2,StressD,[90,100,110,120,140,160,180,200,220],cmap='Greens')

##Plotting final TL from this study for CONUS (Figure 14)############################################################################################
#cont = ax.contourf(xx2, yy2, TL2,[2,4,6,8,10,12,14,16],norm=TwoSlopeNorm(vmin = 0, vcenter = 10, vmax = 20), cmap='Blues')

##Plotting the difference between approx code and our TL (Figure 16)#################################################################################
cont = ax.contourf(xx2, yy2, TLC2-TL2,[-7,-5,-3,-2,-1,-0.5,0.5,1,2,3],norm=TwoSlopeNorm(vmin = -7, vcenter = 0, vmax = 7),cmap='coolwarm_r')

## plot colorbar###############################################################
plt.rc('axes', labelsize=45)
ax.tick_params(axis='both',which='major',width=2,length=10)
for c in cont.collections:
      c.set_clip_path(patch)
cbar = plt.colorbar(cont,fraction=0.02, pad=0.02)
cbar.ax.tick_params(axis='both',which='major',width=2,length=10)
################uncomment 1 of the cbar.set_label lines########################
###############################################################################
##Plotting modal magnitude difference at different spectral periods for CONUS (Figure 8a)############################################################
#cbar.set_label("$M_w$$_($$_4$$_s$$_)$-$M_w$$_($$_2$$_s$$_)$",rotation=270,labelpad=60) #50

##Plotting TL difference at different spectral periods for CONUS (Figure 8b)#########################################################################
#cbar.set_label("$T_L$$_($$_4$$_s$$_)$- $T_L$$_($$_2$$_s$$_)$", rotation=270,labelpad=60) #50

##Plotting modal magnitude at spectral period of 2s for CONUS (Figure 9)###########################################################################
#cbar.set_label('Modal Magnitude, $M_w$',rotation=270,labelpad=60)

##Plotting stress drop for CONUS (Figure 12) ########################################################################################################
#cbar.set_label('Stress Drop (bars)',rotation=270,labelpad=50)

##Plotting final TL from this study for CONUS (Figure 14)############################################################################################
#cbar.set_label("$T_L$", rotation=270,labelpad=35)

##Plotting the difference between approx code and our TL (Figure 16)#################################################################################
cbar.set_label("$T_L$$_($$_N$$_E$$_H$$_R$$_P$$_)$$_s$$_i$$_m$$_.$ - $T_L$", rotation=270,labelpad=45)

##############################################################################
##                          Plotting Shapefiles                           ##
##############################################################################
#plotting shapefile, found code snippet at: https://chrishavlin.com/2016/11/16/shapefiles-tutorial/
#Lower 48 states
sff = shapefile.Reader('cb_2018_us_state_5m')
for shape in list(sff.iterShapes()):
    npoints=len(shape.points) # total points
    nparts = len(shape.parts) # total parts
    
    if nparts == 1:
        x_lon = np.zeros((len(shape.points),1))
        y_lat = np.zeros((len(shape.points),1))
        coord = np.zeros((len(shape.points),1))
        for ip in range(len(shape.points)):
            x_lon[ip] = shape.points[ip][0]
            y_lat[ip] = shape.points[ip][1] 
        plt.plot(x_lon,y_lat,color='gray',linewidth=1)
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
            plt.plot(x_lon,y_lat,color='gray',linewidth=1)

#Border of CONUS
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


plt.xlim(xmin-2,xmax)
plt.ylim(ymin,ymax)
plt.rc('axes', labelsize=40)
plt.rc('font',size=40)
csfont = {'fontname':'Times New Roman'}
ax.set_xlabel('Longitude',fontsize=45)
ax.set_ylabel('Latitude',fontsize=45)
plt.show()

#plot clipped points
#t_l_points = ax.scatter(xc,yc,s=0.5,color='k')
#t_l_points.set_clip_path(patch)

