import numpy as np
import matplotlib.pyplot as plt
import os

#Set example number
example=1002

path_figures = '/Users/gretacappello/Desktop/PROJECT_2_METIS_TS/PSI_simulation_erika/plots/'
folderr = os.path.join(path_figures, f"test_{example}")
cartella = folderr + "/"
npz_path = os.path.join(cartella, "brightness_results.npz")

data = np.load(npz_path, allow_pickle=True)

angles = data["angles"]
c_values = data["c_values"]
delta_B_all = data["delta_B_all"]
delta_B_min = data["delta_B_min"]
delta_B_max = data["delta_B_max"]
a_Bmin = data["a_Bmin"]
c_Bmin = data["c_Bmin"]
a_Bmax = data["a_Bmax"]
c_Bmax = data["c_Bmax"]
sigma_min_limit = data["sigma_min_limit"].item() if "sigma_min_limit" in data else min(delta_B_min)

# --- Prepare figure ---
fig, ax = plt.subplots(figsize=(9,5), dpi=100)

# --- Scatter all delta_B values for each c ---
for idx, c_v in enumerate(c_values):
    delta_B = delta_B_all[idx]
    plt.scatter(angles, c_v*np.ones_like(angles), c=delta_B, s=3, cmap='Oranges_r', vmin=sigma_min_limit, vmax=10*sigma_min_limit)

# --- Highlight regions with low sigma_B ---
B_array = np.array(delta_B_min)
g_array = np.array(a_Bmin)
c_array = np.array(c_Bmin)
n = np.where(B_array < 2 * sigma_min_limit)

res_gamma = a_Bmin[delta_B_min.tolist().index(min(delta_B_min))]
res_c = c_Bmin[delta_B_min.tolist().index(min(delta_B_min))]

# --- Customize plot ---
plt.yscale('log') 
plt.colorbar(label='$\\sigma_B$ (MSB)')

# Vertical lines: gamma
plt.axvline(res_gamma, color='dimgray', linestyle='--', zorder=1)
plt.axvline(min(g_array[n]), color='darkgray', linestyle='--', zorder=1)
plt.axvline(max(g_array[n]), color='darkgray', linestyle='--', zorder=1)

# Horizontal lines: C
plt.axhline(res_c, color='dimgray', linestyle='--', zorder=1)
plt.axhline(min(c_array[n]), color='darkgray', linestyle='--', zorder=1)
plt.axhline(max(c_array[n]), color='darkgray', linestyle='--', zorder=1)

# Scatter points
plt.scatter(g_array[n], c_array[n], marker='s', color='k', s=8, label='$\\sigma_B$ < $2\\cdot$min $\\sigma_B$', zorder=2)
plt.scatter(res_gamma, res_c, marker='o', color='white', edgecolors='tab:orange', label='min $\\sigma_B$', zorder=2)

# Labels and legend
plt.xlabel('$\\gamma$ (°)')
plt.ylabel('$C$ (MSB m$^2$)')
plt.legend()

# Axes style
ax.xaxis.set_tick_params(direction='in', which='both')
ax.yaxis.set_tick_params(direction='in', which='both')

# Optional title with uncertainty bounds
# plt.title(
#     f"$\\gamma = {int(res_gamma)}^{{+{max(g_array[n]) - res_gamma:.2f}}}_{{-{res_gamma - min(g_array[n]):.2f}}}$ "
#     f"and C = {res_c:.2e}^{{+{max(c_array[n]) - res_c:.2e}}}_{{-{res_c - min(c_array[n]):.2e}}} \\, \\mathrm{{MSB}} \\, \\mathrm{{m}}^2$"
# )

# --- Save or show ---

np.savez("brightness_fit_data.npz",
         x_coords=x_coords,
         y_coords=y_coords,
         elong=elong,
         ray_cal=ray_cal,
         res_gamma=res_gamma,
         res_c=res_c,
         alpha=alpha,
         d=d,
         g_array=g_array,
         n=n,
         m_seq_data=m_seq[img_nm].data)
plot_path = os.path.join(cartella, "BRIGHTNESS_trend_REPLOT.png")
plt.savefig(plot_path, dpi=300)
plt.show()
