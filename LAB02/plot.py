import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

t, x, y, dt, E, L = np.loadtxt('out_RK4.dat', unpack=True)

# Trajektoria
plt.plot(x, y, label='trajectory')
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Trajectory')
plt.grid(True)
plt.ticklabel_format(style='plain')
plt.savefig('RK4_trajectory.png', dpi=300, bbox_inches='tight')
plt.close()

# Krok czasowy
plt.plot(t, dt, label='time_step')
plt.legend()
plt.xlabel('t')
plt.ylabel('dt(t)')
plt.title('time_step')
plt.grid(True)
plt.ticklabel_format(style='plain')
plt.savefig('RK4_time_step.png', dpi=300, bbox_inches='tight')
plt.close()

# Energia
plt.plot(t, E, label='total Energy')
plt.legend()
plt.xlabel('t')
plt.ylabel('E(t)')
plt.title('total Energy')
plt.grid(True)
plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter('%.8f'))
plt.savefig('RK4_total_energy.png', dpi=300, bbox_inches='tight')
plt.close()

# Moment pędu
plt.plot(t, L, label='momentum')
plt.legend()
plt.xlabel('t')
plt.ylabel('L(t)')
plt.title('momentum')
plt.grid(True)
plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter('%.8f'))
plt.savefig('RK4_momentum.png', dpi=300, bbox_inches='tight')
plt.close()
