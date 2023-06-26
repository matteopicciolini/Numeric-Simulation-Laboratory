import matplotlib.pyplot as plt
import glob
import imageio
import numpy as np
import os

# Percorso dei file dati
file_pattern = "../Data/10.1_city_coord_true_american_0_*.dat"

# Lista per memorizzare i frame dell'animazione
frames = []

# Carica i dati dai file e crea i frame dell'animazione
for file_name in sorted(glob.glob(file_pattern)):
    
    data = np.loadtxt(file_name)
    # Crea il grafico
    fig, ax = plt.subplots()
    ax.plot(data[:, 0], data[:, 1], marker = "o", linestyle = '-', lw = 1, markersize = 4)
    ax.set_title(f"Generation {file_name[-7:-4]}")
    ax.plot(data[:, 0][0], data[:, 1][0], marker = "X", c = 'darkorange', markersize = 10)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.grid(True)
    plt.savefig("temp_frame.png")
    frames.append(imageio.imread("temp_frame.png"))
    plt.clf()

# Salva i frame come animazione GIF
imageio.mimsave("animation.gif", frames, duration=0.05)
os.remove("temp_frame.png")
