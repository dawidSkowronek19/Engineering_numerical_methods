import numpy as np
import matplotlib.pyplot as plt
import os


os.makedirs("graphs", exist_ok=True)

def get_grid_data(filename):

    data = np.genfromtxt(filename)

    data = data[~np.isnan(data).any(axis=1)]
    
    x = np.unique(data[:, 0])
    y = np.unique(data[:, 1])

    val = data[:, 2].reshape((len(x), len(y)))
    return x, y, val

def plot_heatmap(filename, title, output_name, cmap='viridis'):

    x, y, val = get_grid_data(filename)
    val_t = val.T
    
    plt.figure(figsize=(10, 4))
    plt.imshow(val_t, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()], 
               cmap=cmap, aspect='auto', interpolation='nearest')
    plt.colorbar(label='Wartość')
    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(f"graphs/{output_name}")
    plt.close()

def plot_error(filename, Q_val):

    if os.path.exists(filename):
        data = np.loadtxt(filename)
        plt.figure(figsize=(8, 5))
        plt.plot(data[:, 0], data[:, 1])
        plt.title(f'Zbieżność dla Q={Q_val}')
        plt.xlabel('Iteracja')
        plt.ylabel('Błąd (Gamma)')
        plt.grid(True)
        plt.savefig(f"graphs/error_Q_{Q_val}.png")
        plt.close()

def plot_streamlines(q):


    x_psi, y_psi, psi = get_grid_data(f"out_psi_Q={q}.dat")

    x_u, y_u, u = get_grid_data(f"out_u_Q={q}.dat")
    _, _, v = get_grid_data(f"out_v_Q={q}.dat")

    u_t = u.T
    v_t = v.T
    psi_t = psi.T
    
    plt.figure(figsize=(10, 4))

    plt.imshow(psi_t, origin='lower', extent=[x_psi.min(), x_psi.max(), y_psi.min(), y_psi.max()], 
               cmap='viridis', aspect='auto', alpha=0.3)

    plt.streamplot(x_u, y_u, u_t, v_t, color='white', density=1.5, linewidth=1)
    
    plt.title(rf"Linie pola prędkości (Q={q})")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(f"graphs/streamlines_Q={q}.png")
    plt.close()

Q_values = [-1000, -4000, 4000]

for q in Q_values:
    plot_streamlines(q)

    plot_heatmap(f"out_zeta_Q={q}.dat", rf"Wirowość $\zeta$ (Q={q})", f"zeta_Q={q}.png", cmap='RdBu_r')

    plot_error(f"out_error_Q={q}.dat", q)
