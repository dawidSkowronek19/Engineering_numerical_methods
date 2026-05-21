import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import os
import re


os.makedirs("graphs", exist_ok=True)

def safe_read(filename):
    data = np.genfromtxt(filename)
    return data[~np.isnan(data).any(axis=1)]

def plot_diagnostics():
    # Wykres średniego położenia X
    if os.path.exists("out_averageX.txt"):
        avg_x = np.loadtxt("out_averageX.txt")
        plt.figure(figsize=(8, 5))
        plt.plot(avg_x[:, 0], avg_x[:, 1], color='blue')
        plt.title("Średnie położenie X w czasie")
        plt.xlabel("Czas t [s]")
        plt.ylabel("X_avg")
        plt.grid(True)
        plt.savefig("graphs/srednie_polozenie_X.png")
        plt.close()


def create_gif():
    files = glob.glob('./u/t=*.txt')

    def get_t_val(f):
        match = re.search(r't=(\d+)\.txt', f)
        return int(match.group(1)) if match else -1
    files.sort(key=get_t_val)

    files = files[::10]

    first_data = safe_read(files[0])
    x = np.unique(first_data[:, 0])
    y = np.unique(first_data[:, 1])
    nx, ny = len(x), len(y)
    u_max = np.max(first_data[:, 2])

    fig, ax = plt.subplots(figsize=(10, 4))
    im = ax.imshow(np.zeros((ny, nx)), origin='lower', 
                   extent=[x.min(), x.max(), y.min(), y.max()],
                   cmap='jet', vmin=0, vmax=u_max, aspect='auto')
    
    fig.colorbar(im, ax=ax, label='Stężenie u')
    title = ax.set_title("Rozkład u")
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    def update(frame_file):
        data = safe_read(frame_file)
        val = data[:, 2].reshape((nx, ny)).T
        im.set_data(val)
        
        it_num = get_t_val(frame_file)
        title.set_text(rf"Rozkład stężenia $u$ (iteracja: {it_num})")
        return [im, title]

    ani = animation.FuncAnimation(fig, update, frames=files, blit=True)
    

    ani.save("graphs/animacja_rozkladu.gif", fps=10)
    plt.close()

if __name__ == "__main__":
    plot_diagnostics()
    create_gif()