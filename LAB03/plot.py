import numpy as np
import matplotlib.pyplot as plt


data_pot = np.genfromtxt('map_of_potential.txt')

data_pot = data_pot[~np.isnan(data_pot).any(axis=1)]

x = data_pot[:, 0]
y = data_pot[:, 1]
V = data_pot[:, 2]


n_x = int(x.max() / 0.1)
n_y = int(y.max() / 0.1)
X = x.reshape((n_x + 1, n_y + 1))
Y = y.reshape((n_x + 1, n_y + 1))
V_grid = V.reshape((n_x + 1, n_y + 1))


data_pos = np.loadtxt('position.txt')
t_pos = data_pos[:, 0]
x_pos = data_pos[:, 1]
y_pos = data_pos[:, 2]


plt.figure(figsize=(10, 8))


contour = plt.contourf(X, Y, V_grid, levels=50, cmap='viridis')
plt.colorbar(contour, label='Potencjał V')


plt.plot(x_pos, y_pos, color='red', linewidth=2, label='Trajektoria cząstki')
plt.scatter(x_pos[0], y_pos[0], color='white', marker='o', label='Start') # Punkt startowy

plt.title('Mapa potencjału i trajektoria cząstki')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(alpha=0.3)

plt.savefig("./output.png")