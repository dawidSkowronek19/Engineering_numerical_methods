import numpy as np 
import matplotlib.pyplot as plt
t, u, z, mi=np.loadtxt('out_picard.dat', unpack=True)

plt.plot(t, u, label='u(t)')
plt.plot(t,z,label='z(t)')
plt.legend()
plt.xlabel('t')
plt.ylabel(' ')
plt.title('Przebieg Choroby')
plt.grid(True)
plt.savefig('Picard.png', dpi=300, bbox_inches='tight')

plt.figure()
t, u, z, mi=np.loadtxt('out_newton.dat', unpack=True)

plt.plot(t, u, label='u(t)')
plt.plot(t,z,label='z(t)')
plt.legend()
plt.xlabel('t')
plt.ylabel(' ')
plt.title('Przebieg Choroby')
plt.grid(True)
plt.savefig('Newton.png', dpi=300, bbox_inches='tight')

plt.figure()
t,u,z,mi=np.loadtxt('out_RK2.dat', unpack=True)
plt.plot(t, u, label='u(t)')
plt.plot(t,z,label='z(t)')
plt.legend()
plt.xlabel('t')
plt.ylabel(' ')
plt.title('Przebieg Choroby')
plt.grid(True)
plt.savefig('RK2.png', dpi=300, bbox_inches='tight')

