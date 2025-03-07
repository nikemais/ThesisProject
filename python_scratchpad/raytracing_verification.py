import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import re

path = r'C:\Users\nike\Documents\ThesisProject\python_scratchpad\data_raytracing\lowfidelity'



with open(path + '\\a_1.txt', 'r') as f:
    n = sum(1 for line in f)
my_c = np.zeros((n, 6))
rt_c = np.zeros((n, 6))
angles = np.zeros((n, 2))

with open(path + '\\c_1.txt', 'r') as f:
    for i in range(0, n):
        my_c[i, :] = f.readline().split('\t')
with open(path + '\\a_1.txt', 'r') as f:
    for i in range(0, n):
        angles[i, :] = f.readline().split('\t')
with open(path + '\\table_v2.txt', 'r') as f:
    for i in range(0, 2701):
        numbers = re.split(r'\s{1,}', f.readline().strip())
        rt_c[i, :] = numbers

delta_beta = round(abs(angles[0, 1] - angles[1, 1])*180/np.pi*100)/100
n_beta = int(360/delta_beta + 1)
n_alpha = n/n_beta
delta_alpha = round(abs(angles[0, 0] - angles[n_beta+1, 0])*180/np.pi*100)/100
alpha = np.arange(-90, 90+delta_alpha, delta_alpha)
beta = np.arange(-180, 180+delta_beta, delta_beta)

# extra rotation for v4: wind-> body
for i in range(n):
    cos_alpha = np.cos(-angles[i, 0])
    cos_beta = np.cos(angles[i, 1]-np.pi/2)
    sin_alpha = np.sin(-angles[i, 0])
    sin_beta = np.sin(angles[i, 1]-np.pi/2)
    R = np.array([cos_alpha*cos_beta, sin_beta, sin_alpha*cos_beta,
                  -cos_alpha*sin_beta, cos_beta, -sin_alpha*sin_beta,
                  -sin_alpha, 0, cos_alpha]).reshape(3, 3)
    my_c[i, 0:3] = R.dot(my_c[i, 0:3])
    my_c[i, 3:6] = R.dot(my_c[i, 3:6])


error = rt_c - my_c
fig, axs = plt.subplots(1, 3, figsize=(15, 5))
values = [r'$\Delta C_{x}$', r'$\Delta C_{y}$', r'$\Delta C_{z}$', r'$C_{x, ir}$', r'$C_{y, ir}$', r'$C_{z, ir}$']
for i in range(3):
    heatmap = error[:, i].reshape(int(n_alpha), n_beta)
    B, A = np.meshgrid(beta, alpha)
    im = axs[i].pcolormesh(B, A, heatmap, shading='auto', cmap='viridis')
    contour_levels = [-0.05, -0.025, 0, 0.025, 0.05]
    cs = axs[i].contour(B, A, heatmap, levels=6, colors='black', linewidths=0.8)
    axs[i].clabel(cs, inline=True, fontsize=8)
    cbar = fig.colorbar(im, ax=axs[i])
    cbar.set_label(values[i] + r'$[m^2$]')

    axs[i].set_xlabel(r'$\beta$ [deg]')
    axs[i].set_ylabel(r'$\alpha$ [deg]')
    axs[i].set_title(values[i])

plt.tight_layout()
plt.savefig('data_raytracing/lowfidelity/error.png', dpi=200)

fig, axs = plt.subplots(2, 3, figsize=(15, 10))
values = [r'$C_{x, op}$', r'$C_{y, op}$', r'$C_{z, op}$', r'$C_{x, ir}$', r'$C_{y, ir}$', r'$C_{z, ir}$']
for i in range(6):
    heatmap = my_c[:, i].reshape(int(n_alpha), n_beta)
    B, A = np.meshgrid(beta, alpha)
    im = axs[i // 3, i % 3].pcolormesh(B, A, heatmap, shading='auto', cmap='viridis')
    contour_levels = [-0.05, -0.025, 0, 0.025, 0.05]
    cs = axs[i // 3, i % 3].contour(B, A, heatmap, levels=6, colors='black', linewidths=0.8)
    axs[i // 3, i % 3].clabel(cs, inline=True, fontsize=8)
    cbar = fig.colorbar(im, ax=axs[i // 3, i % 3])
    cbar.set_label(values[i] + r'$[m^2$]')

    axs[i // 3, i % 3].set_xlabel(r'$\beta$ [deg]')
    axs[i // 3, i % 3].set_ylabel(r'$\alpha$ [deg]')
    axs[i // 3, i % 3].set_title(values[i])

plt.tight_layout()
plt.savefig('data_raytracing/lowfidelity/ssh_results_v1_cl.png', dpi=200)

# fig, axs = plt.subplots(2, 3, figsize=(15, 10))
# values = [r'$C_{x, op}$', r'$C_{y, op}$', r'$C_{z, op}$', r'$C_{x, ir}$', r'$C_{y, ir}$', r'$C_{z, ir}$']
# for i in range(6):
#     heatmap = rt_c[:, i].reshape(37, 73)
#     B, A = np.meshgrid(beta, alpha)
#     im = axs[i // 3, i % 3].pcolormesh(B, A, heatmap, shading='auto', cmap='viridis')
#
#     # Add colorbar with label
#     cbar = fig.colorbar(im, ax=axs[i // 3, i % 3])
#     cbar.set_label(values[i] + r'$[m^2$]')
#
#     # Set axis labels
#     axs[i // 3, i % 3].set_xlabel(r'$\beta$ [deg]')
#     axs[i // 3, i % 3].set_ylabel(r'$\alpha$ [deg]')
#     axs[i // 3, i % 3].set_title(values[i])
#
# plt.tight_layout()
# plt.savefig('rt_results.png', dpi=200)

fig, axs = plt.subplots(2, 3, figsize=(15, 10))
values = [r'$C_{x, ssh}$', r'$C_{y, ssh}$', r'$C_{z, ssh}$']
values_rt = [r'$C_{x, rt}$', r'$C_{y, rt}$', r'$C_{z, rt}$']
# comparison rt and mine
for i in range(3):
    heatmap = my_c[:, i].reshape(int(n_alpha), n_beta)
    B, A = np.meshgrid(beta, alpha)
    im = axs[0, i % 3].pcolormesh(B, A, heatmap, shading='auto', cmap='viridis')

    contour_levels = [-0.05, -0.025, 0, 0.025, 0.05]
    cs = axs[i // 3, i % 3].contour(B, A, heatmap, levels=6, colors='black', linewidths=0.8)
    axs[i // 3, i % 3].clabel(cs, inline=True, fontsize=8)
    cbar = fig.colorbar(im, ax=axs[0, i % 3])
    cbar.set_label(values[i] + r'$[m^2$]')

    axs[0, i % 3].set_xlabel(r'$\beta$ [deg]')
    axs[0, i % 3].set_ylabel(r'$\alpha$ [deg]')
    axs[0, i % 3].set_title(values[i])

    heatmap = rt_c[:, i].reshape(int(n_alpha), n_beta)
    B, A = np.meshgrid(beta, alpha)
    im = axs[1, i % 3].pcolormesh(B, A, heatmap, shading='auto', cmap='viridis')

    cs = axs[1, i % 3].contour(B, A, heatmap, levels=6, colors='black', linewidths=0.8)
    axs[1, i % 3].clabel(cs, inline=True, fontsize=8)
    cbar = fig.colorbar(im, ax=axs[1, i % 3])
    cbar.set_label(values_rt[i] + r'$[m^2$]')

    axs[1, i % 3].set_xlabel(r'$\beta$ [deg]')
    axs[1, i % 3].set_ylabel(r'$\alpha$ [deg]')
    axs[1, i % 3].set_title(values_rt[i])

plt.tight_layout()
plt.savefig('data_raytracing/lowfidelity/comparison_v1_cl.png', dpi=200)
print(np.max(my_c[:, 0]))
print(np.max(rt_c[:, 0]))

fig, axs = plt.subplots(3, 3, figsize=(15, 15))
values = [r'$C_{x, ssh}$', r'$C_{y, ssh}$', r'$C_{z, ssh}$']
values_rt = [r'$C_{x, rt}$', r'$C_{y, rt}$', r'$C_{z, rt}$']
values_err = [r'$\Delta C_{x, rt}$', r'$\Delta C_{y, rt}$', r'$\Delta C_{z, rt}$']
# comparison rt and mine
for i in range(3):
    heatmap = my_c[:, i].reshape(int(n_alpha), n_beta)
    B, A = np.meshgrid(beta, alpha)
    im = axs[0, i % 3].pcolormesh(B, A, heatmap, shading='auto', cmap='viridis')

    contour_levels = [-0.05, -0.025, 0, 0.025, 0.05]
    cs = axs[i // 3, i % 3].contour(B, A, heatmap, levels=6, colors='black', linewidths=0.8)
    axs[i // 3, i % 3].clabel(cs, inline=True, fontsize=8)
    cbar = fig.colorbar(im, ax=axs[0, i % 3])
    cbar.set_label(values[i] + r'$[m^2$]')

    axs[0, i % 3].set_xlabel(r'$\beta$ [deg]')
    axs[0, i % 3].set_ylabel(r'$\alpha$ [deg]')
    axs[0, i % 3].set_title(values[i])

    heatmap = rt_c[:, i].reshape(int(n_alpha), n_beta)
    B, A = np.meshgrid(beta, alpha)
    im = axs[1, i % 3].pcolormesh(B, A, heatmap, shading='auto', cmap='viridis')

    cs = axs[1, i % 3].contour(B, A, heatmap, levels=6, colors='black', linewidths=0.8)
    axs[1, i % 3].clabel(cs, inline=True, fontsize=8)
    cbar = fig.colorbar(im, ax=axs[1, i % 3])
    cbar.set_label(values_rt[i] + r'$[m^2$]')

    axs[1, i % 3].set_xlabel(r'$\beta$ [deg]')
    axs[1, i % 3].set_ylabel(r'$\alpha$ [deg]')
    axs[1, i % 3].set_title(values_rt[i])

    heatmap = error[:, i].reshape(int(n_alpha), n_beta)
    B, A = np.meshgrid(beta, alpha)
    im = axs[2, i % 3].pcolormesh(B, A, heatmap, shading='auto', cmap='viridis')

    cs = axs[2, i % 3].contour(B, A, heatmap, levels=3, colors='black', linewidths=0.8)
    axs[2, i % 3].clabel(cs, inline=True, fontsize=8)
    cbar = fig.colorbar(im, ax=axs[2, i % 3])
    cbar.set_label(values_err[i] + r'$[m^2$]')

    axs[2, i % 3].set_xlabel(r'$\beta$ [deg]')
    axs[2, i % 3].set_ylabel(r'$\alpha$ [deg]')
    axs[2, i % 3].set_title(values_err[i])

plt.tight_layout()
plt.savefig('data_raytracing/lowfidelity/comparison.png', dpi=200)


# test beta = 0 behaviour
# dt = 0.1
# beta_test = np.arange(0, 90+dt, dt)
# A1 = 2*4
# A2 = 2*2
# Ca = 0.9
# Cd = 0.07
# Cs = 0.03


# def a_y(beta_):
#     return (-A1*np.cos(beta_)*(-(Ca+Cd)-(2/3+2*np.cos(beta_)*np.cos(beta_)))-
#             A2*np.sin(beta_)*(-(Ca+Cd)-(2/3+2*np.sin(beta_)*np.sin(beta_))))
#
# ay = np.zeros(len(beta_test))
#
# for i in range(len(beta_test)):
#     ay[i] = a_y(beta_test[i]/180*np.pi)
# plt.figure()
# plt.plot(beta_test, ay)
# plt.savefig('test.png', dpi=200)
