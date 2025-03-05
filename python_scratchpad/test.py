import numpy as np
import matplotlib.pyplot as plt

path = r'C:\Users\nike\Documents\ThesisProject\cpp_scratchpad\output3\PIP_times'
times_ray = np.zeros((22, 100))
times_wind = np.zeros((22, 100))
n = np.zeros(22)
for i in range(0, 100):
    with open(path+str(i+1)+'.txt', 'r') as f:
        for j in range(0, 22):
            line = f.readline().split(' ')
            n[j] = float(line[0])
            times_ray[j, i] = float(line[1])
            times_wind[j, i] = float(line[2])
# average value
times_ray_avg = np.asarray([np.average(times_ray[i, :]) for i in range(0, 22)])
times_wind_avg = np.asarray([np.average(times_wind[i, :]) for i in range(0, 22)])


plt.figure()
plt.loglog(n, times_ray_avg, label='PIP raycast')
plt.loglog(n, times_wind_avg, label='PIP winding')
plt.xlabel('Number of points [-]')
plt.ylabel('Processing time [ms]')
plt.legend()
plt.grid()
plt.savefig('time_benchmark3.png', dpi=100)



