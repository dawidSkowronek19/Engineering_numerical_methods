import numpy as np
import matplotlib.pyplot as plt
import glob


k_values = [16, 8, 4, 2, 1]


fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

for idx, k in enumerate(k_values):
    filename = f"potential_k={k}.out"
    data = np.genfromtxt(filename)

    data = data[~np.isnan(data).any(axis=1)]
    
    x = np.unique(data[:, 0])
    y = np.unique(data[:, 1])
    V = data[:, 2].reshape((len(x), len(y)))
    
    ax = axes[idx]

    im = ax.imshow(V.T, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()], 
                   cmap='viridis', aspect='auto', interpolation='nearest')
    
    ax.set_title(f'Potencjał (k={k})')
    plt.colorbar(im, ax=ax)

axes[5].axis('off')
plt.tight_layout()
plt.savefig("mapy_potencjalu.png")

plt.figure(figsize=(10, 6))
for k in k_values:
    filename = f"error_k={k}.out"
    data = np.genfromtxt(filename)
    plt.semilogy(data[:, 0], data[:, 1], label=f'k={k}')

plt.title('Zbieżność metody')
plt.xlabel('Numer iteracji')
plt.ylabel('S (log)')
plt.legend()
plt.grid(True, which="both", ls="-", alpha=0.5)
plt.savefig("wykres_bledu.png")


plt.savefig("./zbieznosc.png")