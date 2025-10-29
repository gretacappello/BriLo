#conda activate brilo
#cd <directory_to_brilo_python.py>
#python brilo_3tracks_newgeometry_v2.py

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
import pandas as pd
import matplotlib.colors as mcolors
from astropy.coordinates import SkyCoord
import re
from sunpy.coordinates import HeliographicCarrington

print('done')
#Change the following path to save the results longitude_results.txt and longitude_results_carr.txt 
open("/Users/gretacappello/Desktop/PROJECT_2_METIS_TS/PSI_simulation_erika/plots/LX_tests/longitude_results.txt", "w").close()
open("/Users/gretacappello/Desktop/PROJECT_2_METIS_TS/PSI_simulation_erika/plots/LX_tests/longitude_results_carr.txt", "w").close()
#replace with ray = #ray, e.g. ray=5 if the name of the folder containing the 3 test folders is called ray5. In that case remove loop.
#In this specific code we assume that we tracked 10 rays, so the following command would fit per each ray the three tracking tests that are stored in the ray's folder
#For clarification, you should have a folder for each ray, e.g. named ray5/, and inside ray5/ you should have three folders, one for each test e.g., test1/, test2/, test3/
#Each test represent a track. We do three tracks per ray to account for errors.
for ray in range(1, 11):   
    data_type = 'LX' #Possible options of data type are: 'L3' and 'LX'
    path = f'/Users/gretacappello/Desktop/PROJECT_2_METIS_TS/PSI_simulation_erika/plots/LX_data_product/ray{ray}*.fits' #replace with your directory to the fits file
    path_figures = f'/Users/gretacappello/Desktop/PROJECT_2_METIS_TS/PSI_simulation_erika/plots/LX_tests/ray{ray}/' #replace with your directory in which the tests folders for a specific ray are saved
    
    print("Ray:", ray)
    longitude = []
    latitude  = []
    longitude_carr = []
    latitude_carr  = []
    chosen_med = []
    chosen_low = []
    chosen_up = []
    chosen_med_carr = []
    chosen_low_carr = []
    chosen_up_carr = []    
    med0 = []
    low0 = []
    up0 = []
    chosen_med_b = []
    chosen_low_b = []
    chosen_up_b = []
    chosen_med_b_carr = []
    chosen_low_b_carr = []
    chosen_up_b_carr = []

    for example in range(1, 4):  # range goes from 1 to 3  (1, 4), because we have 3 tests per ray.
        
        print("Current value of example:", example)

        key = 'interpolated_'

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

        print("long test psp: ", header[0]['HGLN_OBS'])


        m_seq = sunpy.map.Map(files, sequence = True)


        img_nm = 0


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
        #plt.show()
        plt.close()


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
        #plt.show()
        plt.close()


        #distance of PSP from the Sun in solar radii

        #d_PSP = header[0]['DSUN_OBS']
        #print('d =', d_PSP/(6.96e8))

        d = header[img_nm]['DSUN_OBS']
        alpha = 0  #reference axis

        print('d =', d/(6.96e8))
        print('alpha =', alpha)

        if alpha < 0:
            alpha = alpha + 360

        print('alpha =', alpha)
        print(example)




        #function to calculate the average total brightness for L3 data and image img in a 7x7 box
        def ray_bright(img):
            mn_brt = []
            
            for i in range(len(x_coords)):
                x = x_coords[i]
                y = y_coords[i]
                
                # Extract the 7x7 region around (x, y)
                sub_img = m_seq[img].data[y-3:y+4, x-3:x+4]
                
                # Check if the region is 7x7 (handle boundaries)
                if sub_img.shape == (7, 7):
                    median_val = np.median(sub_img)
                else:
                    median_val = np.nan  # or some default/fallback value
                
                mn_brt.append(median_val)
            
            return mn_brt



        #calculate total brightness

        ray_cal = ray_bright(img_nm)

        print('measured brightness:', ray_cal)
        print('number of data points:', np.size(ray_cal))




        #read from the pixel coordinates the elongation in arcsec and conversion to degrees

        lon_lat = m_seq[img_nm].pixel_to_world(x_coords*u.pixel,y_coords*u.pixel)
        lon_pix = np.deg2rad(lon_lat.Tx.degree)   # "elongation" in radians
        lat_pix = np.deg2rad(lon_lat.Ty.degree)   # latitude in radians

        # Define Sun's reference coordinates (0,0) in radians
        lat1 = 0.0
        lon1 = 0.0

        # Delta longitude
        delta_lon = lon_pix - lon1

        # Central angle (spherical law of cosines) # elong is \varepsilon' in Cappello+ paper
        elong = np.rad2deg(np.arccos(np.sin(lat1) * np.sin(lat_pix) + np.cos(lat1) * np.cos(lat_pix) * np.cos(delta_lon)))

        print('elongation 1 (central angle):', elong)




        #plot measured brightness against the elongation

        plt.figure(dpi=100)
        plt.plot(elong, ray_cal,'-', color = 'tab:blue', alpha = 0.3)
        plt.plot(elong, ray_cal,'o', color = 'tab:blue')
        plt.xlabel('$\\xi$ (°)')
        plt.ylabel('$B$ (MSB)')
        plt.title('Brightness profile')
        plot_path = os.path.join(cartella, "BRIGHTNESS_trend.png")
        plt.savefig(plot_path, dpi=300)
        #plt.show()
        plt.close()
        print(example)




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

        c1 = np.arange(1.e5,1.e6,.1e5)
        c2 = np.arange(.1e6,1.e7,.1e6)
        c3 = np.arange(.1e7,1.e8,.1e7)
        c4 = np.arange(.1e8,1.e9,.1e8)
        c5 = np.arange(.1e9,1.e10,.1e9)
        c6 = np.arange(.1e10,1.e11,.1e10)
        c7 = np.arange(.1e11,1.e12,.1e11)
        c8 = np.arange(.1e12,1.e13,.1e12)


        #c values range would change between L3 and LX dataset because the MSB covers a different range in the two datasets.

        if data_type == 'L3':
                c_values = np.concatenate((np.array(c2), np.array(c3), np.array(c4)))
        if data_type == 'LX':
                c_values = np.concatenate((np.array(c3), np.array(c4), np.array(c5)))


        print(c_values)
        print(np.size(c_values))



        print(example)



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
        #plt.show()
        plt.close()



        #results for gamma and C

        res_gamma = a_Bmin[delta_B_min.index(min(delta_B_min))]
        res_c = c_Bmin[delta_B_min.index(min(delta_B_min))]

        print('Results: ')
        print('sigma = ', min(delta_B_min))
        print('gamma = ', res_gamma , "°")
        print('C = {:.2e}'.format(res_c), 'MSB m^2')
        print(example)



        #same plot with limits set according to sigma_min (sigma_min,10xsigma_min)

        sigma_min_limit = min(delta_B_min)

        step = 0

        delta_B_min = []
        delta_B_max = []
        a_Bmin = []
        c_Bmin = []
        a_Bmax = []
        c_Bmax = []
        delta_B_all = []

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
                
            delta_B_all.append(delta_B)
            delta_B_min.append(min(delta_B))
            a_Bmin.append(a_gloop[np.argmin(delta_B)])
            c_Bmin.append(c_gloop[np.argmin(delta_B)])
            
            delta_B_max.append(max(delta_B))
            a_Bmax.append(a_gloop[np.argmax(delta_B)])
            c_Bmax.append(c_gloop[np.argmax(delta_B)])
            
            plt.scatter(angles, c_v*np.ones(np.size(angles)), c=delta_B, s=3, cmap='Oranges_r')
            plt.clim(sigma_min_limit, 10*sigma_min_limit) #set limits here
            
        B_array = np.array(delta_B_min)
        n = np.where(B_array<2*sigma_min_limit)
        g_array = np.array(a_Bmin)
        c_array = np.array(c_Bmin)

        np.savez(cartella + "brightness_results.npz",
                angles=np.array(angles),
                c_values=np.array(c_values),
                delta_B_all=np.array(delta_B_all, dtype=object),
                delta_B_min=np.array(delta_B_min),
                delta_B_max=np.array(delta_B_max),
                a_Bmin=np.array(a_Bmin),
                c_Bmin=np.array(c_Bmin),
                a_Bmax=np.array(a_Bmax),
                c_Bmax=np.array(c_Bmax),
                sigma_min_limit=sigma_min_limit)

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


        plot_path = os.path.join(cartella, key + "_BRIGHTNESS_trend.png")
        plt.savefig(plot_path, dpi=300)
        plt.close()



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

        #plot best fit curve

        elong_cal = np.linspace(elong[0],elong[-1],100)
        theo_cal = theo_bright(d,np.deg2rad(elong_cal),np.deg2rad(res_gamma),np.deg2rad(alpha),res_c)

        theo_cal_min = theo_bright(d,np.deg2rad(elong_cal),np.deg2rad(min(g_array[n])),np.deg2rad(alpha),res_c)
        theo_cal_max = theo_bright(d,np.deg2rad(elong_cal),np.deg2rad(max(g_array[n])),np.deg2rad(alpha),res_c)

        #results:
        gamma1 = np.deg2rad(res_gamma)        # scalar
        eps = np.deg2rad(lon_lat.Tx.degree)   # array
        beta = np.deg2rad(lon_lat.Ty.degree)  # array
        r1 = header[0]['DSUN_OBS']            # scalar
        elong_rad = np.deg2rad(elong)
        r2 = r1 * np.sin(elong_rad)/np.sin(gamma1+elong_rad)       # array
        print("start check")
        print("eps1 = ", elong)
        print("eps = ", np.rad2deg(eps))

        # --- Shape checks ---
        shapes = [eps.shape, beta.shape, r2.shape]
        if not all(s == shapes[0] for s in shapes):
            raise ValueError(f"Shape mismatch: eps {eps.shape}, beta {beta.shape}, r2 {r2.shape}")

        # --- Computations --- d_PSP_P is d_PSP_P' in Cappello et al paper
        d_PSP_P = r2 * np.sin(gamma1) / np.sin(elong_rad)
        print("d_PSP_P': ", d_PSP_P)

        # Ensure d_PSP_P has same shape as beta and r2
        if d_PSP_P.shape != beta.shape:
            raise ValueError(f"d_PSP_P shape {d_PSP_P.shape} does not match beta shape {beta.shape}")

        #delta2 = np.arcsin((d_PSP_P * np.sin(beta)) / r2)
        delta2 = np.arctan2(d_PSP_P * np.sin(beta), np.sqrt(r2**2 - (d_PSP_P * np.sin(beta))**2))
        delta2_v2 = np.pi - delta2

        # Final gamma
        #gamma_proj = -eps + np.arcsin((r1 * np.sin(eps)) / (r2 * np.cos(delta2)))
        #print("check arctg: ", np.arctan2(r1 * np.sin(eps), np.sqrt((r2 * np.cos(delta2))**2 - (r1 * np.sin(eps))**2))  )
        print("check eps: ", eps)
        #gamma_proj = -eps + np.arctan2(r1 * np.sin(eps), np.sqrt((r2 * np.cos(delta2))**2 - (r1 * np.sin(eps))**2))
        
        # Given arrays eps, beta, delta2, elong_rad, and scalars r1, r2 already shaped/checked
        A = d_PSP_P*np.cos(beta)*np.sin(eps) #r1 * np.sin(eps)
        B = r2 * np.cos(delta2)

        # Numerical safety: clamp small negatives to zero due to roundoff
        denom_sq = B**2 - A**2
        denom_sq = np.clip(denom_sq, 0.0, None)
        root = np.sqrt(denom_sq)

        # Base interior angle in [0, pi/2] via atan2 keeps the quadrant based on cos>=0 branch
        theta0 = np.arctan2(A, root)  # one valid solution
        theta1 = np.pi - theta0        # the other valid solution

        # Your projected gamma (apply the -eps shift)
        gamma_a = theta0
        gamma_b = theta1

        # Normalize to (-pi, pi] if helpful
        #def wrap_pi(x):
        #    return (x + np.pi) % (2*np.pi) - np.pi
        print("*******************")

        print("gamma_a:", np.rad2deg(gamma_a))
        print("gamma_b:", np.rad2deg(gamma_b))

        def wrap_0_2pi(x):
            return x % (2 * np.pi)
        gamma_a = wrap_0_2pi(gamma_a)
        gamma_b = wrap_0_2pi(gamma_b)
        print("gamma_a wrapped:", np.rad2deg(gamma_a))
        print("gamma_b wrapped:", np.rad2deg(gamma_b))
        print("*******************")

        
        print("gamma solution 1 (rad): ", gamma_a)
        print("gamma solution 2 (rad): ", gamma_b)
        print("delta2 solution 1 (rad): ", delta2)
        print("delta2 solution 2 (rad): ", delta2_v2)

        print("gamma solution 1 (deg): ", np.rad2deg(gamma_a))
        print("gamma solution 2 (deg): ", np.rad2deg(gamma_b))
        print("delta2 (deg): ", np.rad2deg(delta2))
        print("delta2 solution 2 (deg): ", np.rad2deg(delta2_v2))
                        
        import numpy as np
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        from sunpy.coordinates import HeliographicStonyhurst, HeliocentricInertial

        v=0
        lat_rayhgs = []
        lon_rayhgs = []
        lat_raycarr= []
        lon_raycarr = []
        for gamma_sol in [gamma_a, gamma_b]:
            # PSP position (r_psp) and velocity (v_psp) in HCI (from header or ephemeris)
            r_psp = np.array([header[0]["HCIX_OBS"], header[0]["HCIY_OBS"], header[0]["HCIZ_OBS"]]) * u.m
            v_psp = np.array([header[0]["HCIX_VOB"],
                                header[0]["HCIY_VOB"],
                                header[0]["HCIZ_VOB"]]) * (u.m / u.s)

            # Build basis
            x_hat = r_psp.value / np.linalg.norm(r_psp.value)   # radial
            z_hat = np.cross(r_psp.value, v_psp.value)
            z_hat /= np.linalg.norm(z_hat)                      # orbital angular momentum
            y_hat = np.cross(z_hat, x_hat)                      # completes triad

            # Now build unit vectors for the rays
            u_vecs = ((np.cos(delta2)*np.cos(gamma_sol))[:,None] * x_hat +
                        (np.cos(delta2)*np.sin(gamma_sol))[:,None] * y_hat +
                        (np.sin(delta2))[:,None] * z_hat)

            coords_hci = SkyCoord(x=u_vecs[:,0]*u.one,
                                    y=u_vecs[:,1]*u.one,
                                    z=u_vecs[:,2]*u.one,
                                    representation_type="cartesian",
                                    frame=HeliocentricInertial,
                                    obstime=header[0]["DATE-OBS"])

            coords_hgs = coords_hci.transform_to(HeliographicStonyhurst)
            coords_carr = coords_hci.transform_to(HeliographicCarrington(obstime=coords_hci.obstime, observer="sun"))

            psp_hgs = SkyCoord(x=r_psp[0], y=r_psp[1], z=r_psp[2],
                                representation_type="cartesian",
                                frame=HeliocentricInertial,
                                obstime=header[0]["DATE-OBS"]).transform_to(HeliographicStonyhurst)

            print(f"Ray HGS longitudes v{v+1}:", coords_hgs.lon.deg)
            print(f"Ray HGS latitude v{v+1}:", coords_hgs.lat.deg)
            lon_rayhgs.append(coords_hgs.lon.deg)
            lat_rayhgs.append(coords_hgs.lat.deg)
            lon_raycarr.append(coords_carr.lon.deg)
            lat_raycarr.append(coords_carr.lat.deg)

            v=v+1

        print("PSP HGS longitude:", psp_hgs.lon.deg)
        print("PSP HGS latitude:", psp_hgs.lat.deg)

        print("HGLN_OBS: ", header[0]["HGLN_OBS"])
        print("HGLT_OBS: ", header[0]["HGLT_OBS"])
        
        print("beta: ", beta)


        def get_median_and_errors(arr):
            median = np.median(arr)
            lower_err = median - np.percentile(arr, 16)
            upper_err = np.percentile(arr, 84) - median
            return median, lower_err, upper_err

        # Compute medians and errors
        med0, low0, up0 = get_median_and_errors(lon_rayhgs[0])
        med1, low1, up1 = get_median_and_errors(lon_rayhgs[1])

        med0b, low0b, up0b = get_median_and_errors(lat_rayhgs[0])
        med1b, low1b, up1b = get_median_and_errors(lat_rayhgs[1])

        med0_carr, low0_carr, up0_carr = get_median_and_errors(lon_raycarr[0])
        med1_carr, low1_carr, up1_carr = get_median_and_errors(lon_raycarr[1])

        med0b_carr, low0b_carr, up0b_carr = get_median_and_errors(lat_raycarr[0])
        med1b_carr, low1b_carr, up1b_carr = get_median_and_errors(lat_raycarr[1])

        print("gamma1: ", gamma1)
        print("gamma_a: ", gamma_a)
        # Define "error range size" as total span
        range0 = np.max(np.abs(np.rad2deg(gamma1) - np.rad2deg(gamma_a)))
        range1 = np.max(np.abs(np.rad2deg(gamma1) - np.rad2deg(gamma_b)))
        #lat_pred_0 = header[0]["HGLT_OBS"] + np.rad2deg(delta2[0])
        #lat_real_0 = np.rad2deg(lat_rayhgs[0])
        #range0 = np.max(np.abs(lat_pred_0 - lat_real_0))
        
        #lat_pred_1 = header[0]["HGLT_OBS"] + np.rad2deg(delta2[1])
        #lat_real_1 = np.rad2deg(lat_rayhgs[1])
        #range1 = np.max(np.abs(lat_pred_1 - lat_real_1))


        #range0 = np.max(np.abs((np.abs(np.rad2deg(delta2[0])) + np.abs(header[0]["HGLT_OBS"]))-np.abs(np.rad2deg(lat_rayhgs[0]))))
        #range1 = np.max(np.abs((np.abs(np.rad2deg(delta2[0])) + np.abs(header[0]["HGLT_OBS"]))-np.abs(np.rad2deg(lat_rayhgs[1]))))
        
        #np.max(np.abs(np.rad2deg(delta2) + np.abs(np.rad2deg(lat_rayhgs[0]))))
        #range1 = #np.max(np.abs(np.rad2deg(delta2) + np.abs(np.rad2deg(lat_rayhgs[1]))))
        
        # Choose solution
        if range0 <= range1:
            longitude.append(lon_rayhgs[0])
            latitude.append(lat_rayhgs[0])
            longitude_carr.append(lon_raycarr[0])
            latitude_carr.append(lat_raycarr[0])
            chosen_med.append(med0)
            chosen_low.append(low0)
            chosen_up.append(up0)
            chosen_med_b.append(med0b)
            chosen_low_b.append(low0b)
            chosen_up_b.append(up0b)
            chosen_med_carr.append(med0_carr)
            chosen_low_carr.append(low0_carr)
            chosen_up_carr.append(up0_carr)
            chosen_med_b_carr.append(med0b_carr)
            chosen_low_b_carr.append(low0b_carr)
            chosen_up_b_carr.append(up0b_carr)
        else:
            longitude.append(lon_rayhgs[1])
            latitude.append(lat_rayhgs[1])
            longitude_carr.append(lon_raycarr[1])
            latitude_carr.append(lat_raycarr[1])
            chosen_med.append(med1)
            chosen_low.append(low1)
            chosen_up.append(up1)
            chosen_med_b.append(med1b)
            chosen_low_b.append(low1b)
            chosen_up_b.append(up1b)
            chosen_med_carr.append(med1_carr)
            chosen_low_carr.append(low1_carr)
            chosen_up_carr.append(up1_carr)
            chosen_med_b_carr.append(med1b_carr)
            chosen_low_b_carr.append(low1b_carr)
            chosen_up_b_carr.append(up1b_carr)

        #print("Chosen longitude median:", chosen_med)   
        #print("Longitude errors: -{:.3f}, +{:.3f}".format(chosen_low, chosen_up))

        np.savez(cartella +"profile_fit_data.npz",
            x_coords=x_coords,
            y_coords=y_coords,
            beta = np.rad2deg(beta),
            eps1=elong,
            eps = np.rad2deg(eps),
            ray_cal=ray_cal,
            error_bars = error_bars,
            elong_cal = elong_cal,
            theo_cal = theo_cal,
            delta2 = np.rad2deg(delta2),
            gamma1 = np.rad2deg(gamma1),
            gamma_a = np.rad2deg(gamma_a), 
            gamma_b = np.rad2deg(gamma_b), 
            gamma_a_hgs = lon_rayhgs[0],
            delta2_a_hgs = lat_rayhgs[0],
            gamma_b_hgs = lon_rayhgs[1],
            delta2_b_hgs = lat_rayhgs[1],
            gamma_a_carr = lon_raycarr[0],
            delta2_a_carr = lat_raycarr[0],
            gamma_b_carr = lon_raycarr[1],
            delta2_b_carr = lat_raycarr[1],
            r1 = r1,
            r2=r2,
            d_PSP_P = d_PSP_P
        )
        
        df = pd.DataFrame({
            "x_coords": x_coords,
            "y_coords": y_coords,
            "beta_deg": np.rad2deg(beta),
            "eps_deg": np.rad2deg(eps),
            "eps1_deg": elong,
            "ray_cal": ray_cal,
            "error_bars": error_bars,
            "elong_cal": elong_cal,
            "theo_cal": theo_cal,
            "r1": r1,
            "r2": r2,
            "d_PSP_P": d_PSP_P,
            "delta2_deg": np.rad2deg(delta2),
            "delta2_deg_sol2": np.rad2deg(delta2_v2),
            "gamma1_deg": np.rad2deg(gamma1), 
            "gamma_a": np.rad2deg(gamma_a), 
            "gamma_b": np.rad2deg(gamma_b), 
            "gamma_a_hgs_deg": lon_rayhgs[0],
            "gamma_b_hgs_deg": lon_rayhgs[1],
            "delta2_a_hgs_deg": lat_rayhgs[0],
            "delta2_b_hgs_deg": lat_rayhgs[1],
            "gamma_a_carr_deg": lon_raycarr[0],
            "gamma_b_carr_deg": lon_raycarr[1],
            "delta2_a_carr_deg": lat_raycarr[0],
            "delta2_b_carr_deg": lat_raycarr[1], 
            "psp_hgs_lon": psp_hgs.lon.deg, 
            "psp_hgs_lat": psp_hgs.lat.deg, 
            "check_long_PSP": header[0]["HGLN_OBS"], 
            "check_lat_PSP": header[0]["HGLT_OBS"]
        })

        # Save CSV in the same folder
        csv_path = os.path.join(cartella, "profile_fit_data_quicklook.csv")
        df.to_csv(csv_path, index=False)

        fig, ax = plt.subplots(figsize=(7.5,5),dpi=100)
        ax.plot(elong_cal, theo_cal, label=(
            r"best-fit curve with" + "\n"
            + fr"$\gamma1={int(res_gamma)}^\circ$" + "\n"
            + fr"$\gamma={int(np.min(coords_hgs.lon.deg))}, {int(np.max(coords_hgs.lon.deg))}^\circ$" + "\n"
            + fr"$\delta_2={int(np.min(coords_hgs.lat.deg))}, {int(np.max(coords_hgs.lat.deg))}^\circ$" + "\n"
            + r"C = {:.1e} MSB m$^2$".format(res_c)
        ), color='k')

        ax.errorbar(elong, ray_cal, error_bars, label = 'measured data', color = 'tab:orange', fmt = 'x', capsize = 3)
        plt.xlabel('$\\xi$ (°)')
        plt.ylabel('B (MSB)')
        plt.legend()

        ax.xaxis.set_tick_params(direction='in', which='both')
        ax.yaxis.set_tick_params(direction='in', which='both')
        plot_path = os.path.join(cartella, key + "_BRIGHTNESS_fit.png")
        plt.savefig(plot_path, dpi=300)
        plt.close()

    if data_type == "LX":
        with open("/Users/gretacappello/Desktop/PROJECT_2_METIS_TS/PSI_simulation_erika/plots/LX_tests/longitude_results.txt", "a") as f:
            f.write(f"Ray {ray} - Data type: {data_type}\n")  # opzionale, per sapere a quale iterazione si riferisce
            for i in [0, 1, 2]:
                f.write(f"Test {i+1}: Chosen longitude median: {chosen_med[i]:.6f} "
                    f"-{chosen_low[i]:.6f}, +{chosen_up[i]:.6f}\n")
                f.write(f"Test {i+1}: Chosen latitude median: {chosen_med_b[i]:.6f} "
                    f"-{chosen_low_b[i]:.6f}, +{chosen_up_b[i]:.6f}\n")
            f.write("\n")  # riga vuota per separare i diversi rays
        with open("/Users/gretacappello/Desktop/PROJECT_2_METIS_TS/PSI_simulation_erika/plots/LX_tests/longitude_results_carr.txt", "a") as f:
            f.write(f"Ray {ray} - Data type: {data_type}\n")  # opzionale, per sapere a quale iterazione si riferisce
            for i in [0, 1, 2]:
                f.write(f"Test {i+1}: Chosen longitude median: {chosen_med_carr[i]:.6f} "
                    f"-{chosen_low_carr[i]:.6f}, +{chosen_up_carr[i]:.6f}\n")
                f.write(f"Test {i+1}: Chosen latitude median: {chosen_med_b_carr[i]:.6f} "
                    f"-{chosen_low_b_carr[i]:.6f}, +{chosen_up_b_carr[i]:.6f}\n")
            f.write("\n")  # riga vuota per separare i diversi rays
        

    else:
        with open("/Users/gretacappello/Desktop/PROJECT_2_METIS_TS/PSI_simulation_erika/plots/L3_tests/longitude_results.txt", "a") as f:
            f.write(f"Ray {ray} - Data type: {data_type}\n")  # opzionale, per sapere a quale iterazione si riferisce
            for i in [0, 1, 2]:
                f.write(f"Test {i+1}: Chosen longitude median: {chosen_med[i]:.6f} "
                    f"-{chosen_low[i]:.6f}, +{chosen_up[i]:.6f}\n")
                f.write(f"Test {i+1}: Chosen latitude median: {chosen_med_b[i]:.6f} "
                    f"-{chosen_low_b[i]:.6f}, +{chosen_up_b[i]:.6f}\n")
            f.write("\n")  # riga vuota per separare i diversi rays
        with open("/Users/gretacappello/Desktop/PROJECT_2_METIS_TS/PSI_simulation_erika/plots/L3_tests/longitude_results_carr.txt", "a") as f:
            f.write(f"Ray {ray} - Data type: {data_type}\n")  # opzionale, per sapere a quale iterazione si riferisce
            for i in [0, 1, 2]:
                f.write(f"Test {i+1}: Chosen longitude median: {chosen_med_carr[i]:.6f} "
                    f"-{chosen_low_carr[i]:.6f}, +{chosen_up_carr[i]:.6f}\n")
                f.write(f"Test {i+1}: Chosen latitude median: {chosen_med_b_carr[i]:.6f} "
                    f"-{chosen_low_b_carr[i]:.6f}, +{chosen_up_b_carr[i]:.6f}\n")
            f.write("\n")  # riga vuota per separare i diversi rays

if data_type=='L3':
    # Read the file
    with open("/Users/gretacappello/Desktop/PROJECT_2_METIS_TS/PSI_simulation_erika/plots/L3_tests/longitude_results.txt", "r") as f:
        text = f.read()

    # Split by rays
    ray_blocks = re.split(r"\n(?=Ray \d+)", text.strip())

    results = []

    for block in ray_blocks:
        lines = block.strip().split("\n")
        ray_id = re.search(r"Ray (\d+)", lines[0]).group(1)

        lon_vals, lat_vals = [], []

        for line in lines[1:]:
            if "longitude median" in line:
                median = float(line.split()[5])  # 6th element is the value
                lon_vals.append(median)
            elif "latitude median" in line:
                median = float(line.split()[5])
                lat_vals.append(median)

        # Average values
        lon_avg = np.mean(lon_vals)
        lat_avg = np.mean(lat_vals)

        # Errors
        lon_err_plus = max(lon_vals) - lon_avg
        lon_err_minus = lon_avg - min(lon_vals)
        lat_err_plus = max(lat_vals) - lat_avg
        lat_err_minus = lat_avg - min(lat_vals)

        results.append([ray_id, lon_avg, lon_err_plus, lon_err_minus,
                                lat_avg, lat_err_plus, lat_err_minus])

    # Put into DataFrame
    df = pd.DataFrame(results, columns=[
        "Ray", "Avg Lon", "+Err Lon", "-Err Lon", "Avg Lat", "+Err Lat", "-Err Lat"
    ])
    print("HGS frame")
    print(df)

    df.to_csv("/Users/gretacappello/Desktop/PROJECT_2_METIS_TS/PSI_simulation_erika/plots/L3_tests/ray_summary.csv", index=False)



else:
    # Read the file
    with open("/Users/gretacappello/Desktop/PROJECT_2_METIS_TS/PSI_simulation_erika/plots/LX_tests/longitude_results.txt", "r") as f:
        text = f.read()

    # Split by rays
    ray_blocks = re.split(r"\n(?=Ray \d+)", text.strip())

    results = []

    for block in ray_blocks:
        lines = block.strip().split("\n")
        ray_id = re.search(r"Ray (\d+)", lines[0]).group(1)

        lon_vals, lat_vals = [], []

        for line in lines[1:]:
            if "longitude median" in line:
                median = float(line.split()[5])  # 6th element is the value
                lon_vals.append(median)
            elif "latitude median" in line:
                median = float(line.split()[5])
                lat_vals.append(median)

        # Average values
        lon_avg = np.mean(lon_vals)
        lat_avg = np.mean(lat_vals)

        # Errors
        lon_err_plus = max(lon_vals) - lon_avg
        lon_err_minus = lon_avg - min(lon_vals)
        lat_err_plus = max(lat_vals) - lat_avg
        lat_err_minus = lat_avg - min(lat_vals)

        results.append([ray_id, lon_avg, lon_err_plus, lon_err_minus,
                                lat_avg, lat_err_plus, lat_err_minus])

    # Put into DataFrame
    df = pd.DataFrame(results, columns=[
        "Ray", "Avg Lon", "+Err Lon", "-Err Lon", "Avg Lat", "+Err Lat", "-Err Lat"
    ])
    print("HGS frame")
    print(df)

    df.to_csv("/Users/gretacappello/Desktop/PROJECT_2_METIS_TS/PSI_simulation_erika/plots/LX_tests/ray_summary.csv", index=False)

if data_type=='L3':
    # Read the file
    with open("/Users/gretacappello/Desktop/PROJECT_2_METIS_TS/PSI_simulation_erika/plots/L3_tests/longitude_results_carr.txt", "r") as f:
        text = f.read()

    # Split by rays
    ray_blocks = re.split(r"\n(?=Ray \d+)", text.strip())

    results = []

    for block in ray_blocks:
        lines = block.strip().split("\n")
        ray_id = re.search(r"Ray (\d+)", lines[0]).group(1)

        lon_vals, lat_vals = [], []

        for line in lines[1:]:
            if "longitude median" in line:
                median = float(line.split()[5])  # 6th element is the value
                lon_vals.append(median)
            elif "latitude median" in line:
                median = float(line.split()[5])
                lat_vals.append(median)

        # Average values
        lon_avg = np.mean(lon_vals)
        lat_avg = np.mean(lat_vals)

        # Errors
        lon_err_plus = max(lon_vals) - lon_avg
        lon_err_minus = lon_avg - min(lon_vals)
        lat_err_plus = max(lat_vals) - lat_avg
        lat_err_minus = lat_avg - min(lat_vals)

        results.append([ray_id, lon_avg, lon_err_plus, lon_err_minus,
                                lat_avg, lat_err_plus, lat_err_minus])

    # Put into DataFrame
    df = pd.DataFrame(results, columns=[
        "Ray", "Avg Lon", "+Err Lon", "-Err Lon", "Avg Lat", "+Err Lat", "-Err Lat"
    ])
    print("Carrington frame")
    print(df)

    df.to_csv("/Users/gretacappello/Desktop/PROJECT_2_METIS_TS/PSI_simulation_erika/plots/L3_tests/ray_summary_carr.csv", index=False)



else:
    # Read the file
    with open("/Users/gretacappello/Desktop/PROJECT_2_METIS_TS/PSI_simulation_erika/plots/LX_tests/longitude_results_carr.txt", "r") as f:
        text = f.read()

    # Split by rays
    ray_blocks = re.split(r"\n(?=Ray \d+)", text.strip())

    results = []

    for block in ray_blocks:
        lines = block.strip().split("\n")
        ray_id = re.search(r"Ray (\d+)", lines[0]).group(1)

        lon_vals, lat_vals = [], []

        for line in lines[1:]:
            if "longitude median" in line:
                median = float(line.split()[5])  # 6th element is the value
                lon_vals.append(median)
            elif "latitude median" in line:
                median = float(line.split()[5])
                lat_vals.append(median)

        # Average values
        lon_avg = np.mean(lon_vals)
        lat_avg = np.mean(lat_vals)

        # Errors
        lon_err_plus = max(lon_vals) - lon_avg
        lon_err_minus = lon_avg - min(lon_vals)
        lat_err_plus = max(lat_vals) - lat_avg
        lat_err_minus = lat_avg - min(lat_vals)

        results.append([ray_id, lon_avg, lon_err_plus, lon_err_minus,
                                lat_avg, lat_err_plus, lat_err_minus])

    # Put into DataFrame
    df = pd.DataFrame(results, columns=[
        "Ray", "Avg Lon", "+Err Lon", "-Err Lon", "Avg Lat", "+Err Lat", "-Err Lat"
    ])
    print("Carrington frame")

    print(df)

    df.to_csv("/Users/gretacappello/Desktop/PROJECT_2_METIS_TS/PSI_simulation_erika/plots/LX_tests/ray_summary_carr.csv", index=False)

    
    
