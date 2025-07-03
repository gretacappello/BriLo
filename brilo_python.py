#conda activate brilo
#cd <directory_to_brilo_python.py>
#python brilo_python.py

import matplotlib.pyplot as plt
import numpy as np

import glob

from astropy.io import fits
import astropy.units as u
from astropy.visualization import ImageNormalize
from astropy.time import Time
from astropy.visualization import time_support

import sunpy.map
import os
import matplotlib.colors as mcolors
from astropy.coordinates import SkyCoord

print('done')


#Set example number, depending on how you called the example during the tracking step.
example =1001

#CHOOSE A KEY
#key = 'fake'
key = 'interpolated_'  #this it the one to use for the L3 data products
#key = 'points_'
#key = 'real_reproj_'
#longitude_rt = 175
#latitude_rt = 21 #'M10'

#Set path to the WISPR fits file used for the track
path = '/Users/gretacappello/sunpy/data/psp_l3_wispr_20181106t153305_v3_1221.fits'
#Set path in which the example folder (e.g., test_{example}) is located
path_figures = '/Users/gretacappello/Desktop/PROJECT_2_METIS_TS/PSI_simulation_erika/plots/'

#LX
#path = '/Users/gretacappello/Downloads/E10_June19Event_LX/20211119/psp_LX_wispr_20211119T150322_V1_1211.fits'

#L3
#path = '/Users/gretacappello/sunpy/data/psp_l3_wispr_20211119t150322_v1_1211.fits'

#L2
#path = '/Users/gretacappello/Downloads/psp_L2_wispr_20211119T150322_V1_1211.fits'

#simulation tb
#path = '/Users/gretacappello/Downloads/wispr_unmasked_tb000320.fts'

#simulation pb
#path = '/Users/gretacappello/Downloads/wispr_unmasked_pb000320.fts'


folderr = os.path.join(path_figures, f"test_{example}")
if not os.path.exists(folderr):
    os.makedirs(folderr)
cartella = folderr + "/"

files = glob.glob(path)
files = sorted(files)

header = []
images = []

for i in range(len(files)):
    header.append([])
    images.append([])
for i in range(len(files)):   
    header[i] = fits.getheader(files[i])
    images[i] = fits.getdata(files[i])


#make a maps' sequence

m_seq = sunpy.map.Map(files, sequence = True)


# In[33]:


#image number

img_nm = 0


# In[34]:


if key == 'real_':
    d_PSP = header[0]['DSUN_OBS']
    long_PSP = header[0]['HGLN_OBS']
    lat_PSP = header[0]['HGLT_OBS']

    print(d_PSP)

    print(long_PSP)

    print(lat_PSP)


# In[35]:


#read coordinates from text file

if key == 'fake':
    coordinates = open(cartella + 'example_'+str(example)+'_line.txt', 'r')
if key == 'real_reproj_':
    coordinates = open(cartella + 'points_reprojectedL3_' + str(example) + '_line.txt')
if key == 'real_':
    coordinates = open(cartella + 'points_L3_' + str(example) + '_line.txt')
if key == 'interpolated_':
    coordinates = open(cartella + 'interpolated_points_L3_' + str(example) + '_line.txt')
    
x_coords = []
y_coords = []

for cd in coordinates:
    cds = cd.split()
    x_coords.append(int(cds[0]))
    y_coords.append(int(cds[1]))

coordinates.close()

print(key)
print('x_coords =',x_coords)
print('y_coords =',y_coords)

print(example)


# In[36]:


from astropy.wcs import WCS
#plot map with coordinates
wcs_wispr_inner = WCS(header[0])
fig = plt.figure()
ax = fig.add_subplot(projection=m_seq[img_nm])
m_seq[img_nm].plot(axes=ax, clip_interval=(2,99)*u.percent)

lon = ax.coords[0]
lat = ax.coords[1]
lat.set_major_formatter('dd')
lon.set_major_formatter('dd')
plot_path = os.path.join(cartella, key+"_image_ray.png")
plt.savefig(plot_path, dpi=300)
plt.show()


# In[37]:


from astropy.wcs import WCS
#plot map with coordinates
wcs_wispr_inner = WCS(header[0])
fig = plt.figure()
ax = fig.add_subplot(projection=m_seq[img_nm])
m_seq[img_nm].plot(axes=ax, clip_interval=(2,99)*u.percent)

ax.plot(x_coords, y_coords, marker = 'x', ls = '', color = 'tab:orange', markersize = 3)
lon = ax.coords[0]
lat = ax.coords[1]
lat.set_major_formatter('dd')
lon.set_major_formatter('dd')
plot_path = os.path.join(cartella, key+"_image_ray_track.png")
plt.savefig(plot_path, dpi=300)
plt.show()


# In[38]:


#distance of PSP from the Sun in solar radii, longitude and latitude in HGS coordinates

d_PSP = header[0]['DSUN_OBS']
#long_PSP = header[0]['HGLN_OBS']
#lat_PSP = header[0]['HGLT_OBS']

print('d =', d_PSP/(6.96e8))
#print('long =', long_PSP)
#print('lat =', lat_PSP)


# In[39]:


#distance from the Sun and heliographic longitude of PSP

d = header[img_nm]['DSUN_OBS']
alpha = 0 #header[img_nm]['HGLN_OBS']

print('d =', d/(6.96e8))
print('alpha =', alpha)

if alpha < 0:
    alpha = alpha + 360

print('alpha =', alpha)
print(example)


# In[40]:


#function to calculate the average total brightness for L3 data and image img in a 7x7 box

def ray_bright(img):
    
    ll = 7*7
    mn_brt = []
    
    for i in range(len(x_coords)):
        
        a = 0
        x = x_coords[i]
        y = y_coords[i]
        
        for j in range(x-3, x+4):
            for k in range(y-3, y+4):
                a = a + m_seq[img].data[k,j]
        
        mn_brt.append(a/ll)
    
    return(mn_brt)


# In[41]:


#calculate aver total brightness

ray_cal = ray_bright(img_nm)

print('measured brightness:', ray_cal)
print('number of data points:', np.size(ray_cal))


# In[42]:


#read from the pixel coordinates the elongation in arcsec and conversion to degrees

lon_lat = m_seq[img_nm].pixel_to_world(x_coords*u.pixel,y_coords*u.pixel)
elong = lon_lat.Tx.degree

#print('helioprojective longitude:', lon_lat.Tx.degree)
#print('helioprojective latitude:', lon_lat.Ty.degree)
print('elongation:', elong)


# In[43]:


#plot measured brightness against the elongation

plt.figure(dpi=100)
plt.plot(elong, ray_cal,'-', color = 'tab:blue', alpha = 0.3)
plt.plot(elong, ray_cal,'o', color = 'tab:blue')
plt.xlabel('$\\xi$ (°)')
plt.ylabel('$B$ (MSB)')
plt.title('Brightness profile')
plot_path = os.path.join(cartella, "BRIGHTNESS_trend.png")
plt.savefig(plot_path, dpi=300)
plt.show()
print(example)


# In[44]:


#function to calculate the theoretical brightness based on Nistico et al. (2020)

def theo_bright(dist,xi,gamma,a,C):
    
    # d = psp_sun distance
    # xi = elongation mesured
    # gamma = longitude of the ray to be calculated
    # alpha = angle between reference axis and psp
    # C = pi * solar brightness * rsun^2 * sigma_t * n * h
    
    cos4 = np.cos(gamma-a+xi)*np.cos(gamma-a+xi)*np.cos(gamma-a+xi)*np.cos(gamma-a+xi)
    B = C*(1-cos4)/(dist*dist*np.sin(xi)*np.sin(xi))
    
    return B


# In[45]:


#define possible angles

angles = np.arange(0, 180) #from long PSP the FoV starts at about +40 till 110


# In[46]:


#define possible C-values

#c1_old = np.arange(1.e5,1.e6,.1e5)
#c2_old = np.arange(.1e6,1.e7,.1e6)
#c3_old = np.arange(.1e7,1.e8,.1e7)
#c4_old = np.arange(.1e8,1.e9,.1e8)
#c5_old = np.arange(.1e9,1.e10,.1e9)
#c6_old = np.arange(.1e10,1.e11,.1e10)
#c7_old = np.arange(.1e11,1.e12,.1e11)
#c8_old = np.arange(.1e12,1.e13,.1e12)


#print(np.arange(1e5,1e6,1e4))
#print(np.arange(1e5,1e7,1e5))
#print(np.arange(1e6,1e8,1e6))
#print(np.arange(1e7,1e9,1e7))
#print(np.arange(1e8,1e10,1e8))


#c_values = np.concatenate((np.array(c3_old), np.array(c4_old), np.array(c5_old), np.array(c6_old), np.array(c7_old), np.array(c8_old)))
#print(np.size(c_values_old))



#TO ADD ALL AT THE END
c1 = np.arange(1.e5,1.e6,.1e5)
c2 = np.arange(.1e6,1.e7,.1e6)
c3 = np.arange(.1e7,1.e8,.1e7)
c4 = np.arange(.1e8,1.e9,.1e8)
c5 = np.arange(.1e9,1.e10,.1e9)
c6 = np.arange(.1e10,1.e11,.1e10)
c7 = np.arange(.1e11,1.e12,.1e11)
c8 = np.arange(.1e12,1.e13,.1e12)

#print(c2)

#L3 data
c_values = np.concatenate((np.array(c2),np.array(c3), np.array(c4), np.array(c5)))



#Simulated data
#c_values = np.concatenate((np.array(c4), np.array(c5),  np.array(c6),  np.array(c7), np.array(c8)))




# In[47]:


print(c_values)
print(np.size(c_values))


# In[48]:


print(example)


# In[49]:


#calculate difference between measured and theoretical brightness and make scatter plot
#x-axis: angles; y-axis: C-values; colorbar: error


step = 0

delta_B_min = []
delta_B_max = []

a_Bmin = []
c_Bmin = []

a_Bmax = []
c_Bmax = []

plt.figure(figsize=(9,5),dpi=100)


for c_v in c_values:
    
    delta_B = []
    a_gloop = []
    c_gloop = []
    
    for g in angles:
        
        step = step + 1
        print("step : ", step, '/', np.size(angles)*np.size(c_values))

        d_B = 0.0
        
        for i in range(len(elong)):
            d_B = d_B + np.abs(ray_cal[i] - theo_bright(d,np.deg2rad(elong[i]),np.deg2rad(g),np.deg2rad(alpha),c_v))
        
        d_B = d_B/np.size(ray_cal)
        delta_B.append(d_B)
        a_gloop.append(g)
        c_gloop.append(c_v)
    
    delta_B_min.append(min(delta_B))
    a_Bmin.append(a_gloop[delta_B.index(min(delta_B))])
    c_Bmin.append(c_gloop[delta_B.index(min(delta_B))])
    
    delta_B_max.append(max(delta_B))
    a_Bmax.append(a_gloop[delta_B.index(max(delta_B))])
    c_Bmax.append(c_gloop[delta_B.index(max(delta_B))])
    
    plt.scatter(angles, c_v*np.ones(np.size(angles)), c=delta_B, s=3, cmap='viridis')


plt.yscale('log')
plt.colorbar(label='$\\sigma_B$ (MSB)')    
plt.scatter(a_Bmin[delta_B_min.index(min(delta_B_min))], c_Bmin[delta_B_min.index(min(delta_B_min))],marker='s', color='red', label = 'min')
plt.scatter(a_Bmax[delta_B_max.index(max(delta_B_max))], c_Bmax[delta_B_max.index(max(delta_B_max))], marker='s', color='black', label = 'max')
plt.xlabel('$\\gamma$ (°)')
plt.ylabel('$C$ (MSB m$^2$)')
plt.legend()
res_gamma = a_Bmin[delta_B_min.index(min(delta_B_min))]
res_c = c_Bmin[delta_B_min.index(min(delta_B_min))] 
plt.title('$\\gamma= $'+str(int(res_gamma))+'° and C = {:.2e} MSB m$^2$'.format(res_c))
plt.show()


# In[50]:


#results for gamma and C

res_gamma = a_Bmin[delta_B_min.index(min(delta_B_min))]
res_c = c_Bmin[delta_B_min.index(min(delta_B_min))]

print('Results: ')
print('sigma = ', min(delta_B_min))
print('gamma = ', res_gamma , "°")
print('C = {:.2e}'.format(res_c), 'MSB m^2')
print(example)


# In[51]:


#same plot with limits set according to sigma_min (sigma_min,10xsigma_min)

sigma_min_limit = min(delta_B_min)

step = 0

delta_B_min = []
delta_B_max = []

a_Bmin = []
c_Bmin = []

a_Bmax = []
c_Bmax = []

fig, ax = plt.subplots(figsize=(9,5),dpi=100)


for c_v in c_values:
    
    delta_B = []
    a_gloop = []
    c_gloop = []
    
    for g in angles:
        
        step = step + 1
        print("step : ", step, '/', np.size(angles)*np.size(c_values))

        d_B = 0.0
        for i in range(len(elong)):
            d_B = d_B + np.abs(ray_cal[i] - theo_bright(d,np.deg2rad(elong[i]),np.deg2rad(g),np.deg2rad(alpha),c_v))
        
        d_B = d_B/np.size(ray_cal)
        delta_B.append(d_B)
        a_gloop.append(g)
        c_gloop.append(c_v)
    
    delta_B_min.append(min(delta_B))
    a_Bmin.append(a_gloop[delta_B.index(min(delta_B))])
    c_Bmin.append(c_gloop[delta_B.index(min(delta_B))])
    
    delta_B_max.append(max(delta_B))
    a_Bmax.append(a_gloop[delta_B.index(max(delta_B))])
    c_Bmax.append(c_gloop[delta_B.index(max(delta_B))])
    
    plt.scatter(angles, c_v*np.ones(np.size(angles)), c=delta_B, s=3, cmap='Oranges_r')
    plt.clim(sigma_min_limit, 10*sigma_min_limit) #set limits here
    
B_array = np.array(delta_B_min)
n = np.where(B_array<2*sigma_min_limit) #
g_array = np.array(a_Bmin)
c_array = np.array(c_Bmin)

plt.yscale('log') 
plt.colorbar(label='$\\sigma_B$ (MSB)')

plt.axvline(res_gamma, color = 'dimgray', linestyle = '--', zorder = 1)
plt.axvline(min(g_array[n]), color = 'darkgray', linestyle = '--', zorder = 1)
plt.axvline(max(g_array[n]), color = 'darkgray', linestyle = '--', zorder = 1)

plt.axhline(res_c, color = 'dimgray', linestyle = '--', zorder = 1)
plt.axhline(min(c_array[n]), color = 'darkgray', linestyle = '--', zorder = 1)
plt.axhline(max(c_array[n]), color = 'darkgray', linestyle = '--', zorder = 1)

plt.scatter(g_array[n], c_array[n], marker='s', color='k',s=8, label = '$\\sigma_B$ < $2\\cdot$min $\\sigma_B$', zorder = 2)
plt.scatter(a_Bmin[delta_B_min.index(min(delta_B_min))], c_Bmin[delta_B_min.index(min(delta_B_min))], marker = 'o', color='white',
            edgecolors = 'tab:orange', label = 'min $\\sigma_B$', zorder = 2)
   
plt.legend()
plt.xlabel('$\\gamma$ (°)')
plt.ylabel('$C$ (MSB m$^2$)')

ax.xaxis.set_tick_params(direction='in', which='both')
ax.yaxis.set_tick_params(direction='in', which='both')

#res_gamma = a_Bmin[delta_B_min.index(min(delta_B_min))]
#res_c = c_Bmin[delta_B_min.index(min(delta_B_min))]
#plt.title('$\\gamma= $'+str(int(res_gamma))+'° $\pm $ and C = {:.2e} MSB m$^2$'.format(res_c))
#plt.title(
#    f"$\\gamma = {int(res_gamma)}^{{+{max(g_array[n]) - res_gamma:.2f}}}_{{-{res_gamma - min(g_array[n]):.2f}}}$ "
#    f"and C = {res_c:.2e}^{{+{max(c_array[n]) - res_c:.2e}}}_{{-{res_c - min(c_array[n]):.2e}}} \, \\mathrm{{MSB}} \, \\mathrm{{m}}^2$"
#)

plot_path = os.path.join(cartella, key + "_BRIGHTNESS_trend.png")
plt.savefig(plot_path, dpi=300)
plt.show()


# In[52]:


#error intervals for gamma and C

print(g_array[n])
print(c_array[n])
print(len(B_array[n]))

print("Results key = ", key)
print('Minimum gamma =', min(g_array[n]))
print('Results gamma =', res_gamma)
print('Maximum gamma =', max(g_array[n]))

print('Minimum C-value =', min(c_array[n]))
print('Result C-value =', res_c)
print('Maximum C-value =', max(c_array[n]))

print(example)


output_file = cartella + "results" + key + ".txt"

# Open file in write mode
with open(output_file, "w") as f:
    # Write results to file
    f.write(f"Results key = {key}\n")
    f.write(f"Minimum gamma = {min(g_array[n])}\n")
    f.write(f"Results gamma = {res_gamma}\n")
    f.write(f"Maximum gamma = {max(g_array[n])}\n\n")

    f.write(f"Minimum C-value = {min(c_array[n])}\n")
    f.write(f"Result C-value = {res_c}\n")
    f.write(f"Maximum C-value = {max(c_array[n])}\n")

print(f"Results saved to {output_file}")


# In[53]:


#calculate errorbars

error_bars = []
    
for i in range(len(x_coords)):
        
        error_bars_cal = []
        
        x = x_coords[i]
        y = y_coords[i]
        
        for j in range(x-3, x+4):
            for k in range(y-3, y+4):
                error_bars_cal.append(m_seq[img_nm].data[k,j])
                
        error_bars.append(np.std(error_bars_cal))


# In[54]:


#plot best fit curve

elong_cal = np.linspace(elong[0],elong[-1],100)
theo_cal = theo_bright(d,np.deg2rad(elong_cal),np.deg2rad(res_gamma),np.deg2rad(alpha),res_c)

theo_cal_min = theo_bright(d,np.deg2rad(elong_cal),np.deg2rad(min(g_array[n])),np.deg2rad(alpha),res_c)
theo_cal_max = theo_bright(d,np.deg2rad(elong_cal),np.deg2rad(max(g_array[n])),np.deg2rad(alpha),res_c)

fig, ax = plt.subplots(figsize=(7.5,5),dpi=100)
ax.plot(elong_cal, theo_cal, label = 'best-fit curve with $\\gamma= $'+str(int(res_gamma))+'° and C = {:.1e} MSB$\,$m$^2$'.format(res_c), color = 'k')


#ax.plot(elong_cal, theo_cal_min, label = 'curves with $\\gamma=120^{\circ+8^\circ}_{-8^\circ}$ and C$=$4.4e+07 MSB$\,$m$^2$', color = 'gray', linestyle = '--')
#ax.plot(elong_cal, theo_cal_max, color = 'gray', linestyle = '--')

#ax.plot(elong, ray_cal,  label = 'measured data', color = 'tab:orange')

ax.errorbar(elong, ray_cal, error_bars, label = 'measured data', color = 'tab:orange', fmt = 'x', capsize = 3)
plt.xlabel('$\\xi$ (°)')
plt.ylabel('B (MSB)')
plt.legend()

#plt.fill_between(elong_cal, theo_cal_min, theo_cal_max, color = 'lightgray', alpha = 0.5)

ax.xaxis.set_tick_params(direction='in', which='both')
ax.yaxis.set_tick_params(direction='in', which='both')
plot_path = os.path.join(cartella, key + "_BRIGHTNESS_fit.png")
plt.savefig(plot_path, dpi=300)
plt.show()


# In[55]:


key


# In[56]:


print(example)


# In[57]:


print(key + "_BRIGHTNESS_trend.png")


