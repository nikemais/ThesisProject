import numpy as np
import matplotlib.pyplot as plt

beta = np.arange(-90, 90, 0.1)/180*np.pi
v = np.array([0, 1])
A1 = A3 = 1
A2 = 2
ca = 1
cd = 0
cs = 0
a = np.zeros((len(beta), 2))
a1 = np.zeros((len(beta), 2))
a2 = np.zeros((len(beta), 2))
a3 = np.zeros((len(beta), 2))
a_test = np.zeros((len(beta), 2))
for i in range(len(beta)):

    n1 = np.array([np.cos(beta[i]), -np.sin(beta[i])])
    n2 = np.array([-np.sin(beta[i]), -np.cos(beta[i])])
    n3 = np.array([-np.cos(beta[i]), np.sin(beta[i])])

    theta_1 = np.heaviside(np.cos(np.pi/2-beta[i]), 0)*np.cos(np.pi/2-beta[i])
    theta_2 = np.heaviside(np.cos(np.pi/2 + beta[i]), 0)*np.cos(np.pi/2+beta[i])

    a[i, :] = (-A1*theta_1*(-v*(ca+cd) + n1*(2/3*cd+cs*2*theta_1)) -
                A2*np.cos(beta[i])*(-v*(ca+cd) + n2*(2/3*cd+cs*2*np.cos(beta[i]))) -
                A3*theta_2*(-v*(ca+cd) + n3*(2/3*cd+cs*2*theta_2)))
    a1[i, :] = -A1*theta_1*(-v*(ca+cd) + n1*(2/3*cd+cs*2*theta_1))
    a2[i, :] = -A2*np.cos(beta[i])*(-v*(ca+cd) + n2*(2/3*cd+cs*2*np.cos(beta[i])))
    a3[i, :] = -A3*theta_2*(-v*(ca+cd) + n3*(2/3*cd+cs*2*theta_2))
    a_test[i, :] = 2*cs*(A1*np.sin(beta[i])**3+A2*np.cos(beta[i])**3)

plt.figure(figsize=(8, 4))
plt.plot(beta*180/np.pi, a[:, 1], label=r'$C_{tot}$')
plt.axvline(np.arctan(A1/A2)*180/np.pi, color='red', linestyle=':')
plt.axvline(-np.arctan(A1/A2)*180/np.pi, color='red', linestyle=':')
plt.plot(beta*180/np.pi, a1[:, 1], label=r'$C_{A1}$')
plt.plot(beta*180/np.pi, a2[:, 1], label=r'$C_{A2}$')
plt.plot(beta*180/np.pi, a3[:, 1], label=r'$C_{A3}$')
plt.title('Pure adsorption (ca=1, cd=cs=0)')
plt.xlabel(r'$\beta$ [deg]')
plt.ylabel(r'$C_y$ [$m^2$]')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('analytical_solution_pureadsorption.png', dpi=200)
pass
