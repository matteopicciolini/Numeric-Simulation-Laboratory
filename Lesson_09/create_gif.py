import matplotlib.pyplot as plt
import glob
import imageio
import numpy as np
import os

# Percorso dei file dati
#file_pattern = "../Data/09.1_square_city_coord_*.dat"
file_pattern = "../Data/09.1_circ_city_coord_*.dat"

# Lista per memorizzare i frame dell'animazione
frames = []

x_sq = [-1, -1, 1, 1, -1]
y_sq = [-1, 1, 1, -1, -1]

theta = np.linspace(0, 2 * np.pi, 100)
x_c = np.cos(theta)
y_c = np.sin(theta)


# Carica i dati dai file e crea i frame dell'animazione
for file_name in sorted(glob.glob(file_pattern)):
    
    data = np.loadtxt(file_name)
    # Crea il grafico
    fig, ax = plt.subplots()
    #ax.plot(x_sq, y_sq, c = 'red')
    ax.plot(x_c, y_c, c = 'red')
    ax.plot(data[:, 0], data[:, 1], marker = "o", linestyle = '-', lw = 1, markersize = 4)
    ax.set_title(f"Generation {file_name[-7:-4]}")
    ax.plot(data[:, 0][0], data[:, 1][0], marker = "X", c = 'darkorange', markersize = 10)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.grid(True)

    ax.set_xticks([-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1])
    ax.set_aspect('equal')
    #ax.set_xlim(-1.3, 1.3)
    #ax.set_ylim(-1.3, 1.3)

    plt.savefig("temp_frame.png")
    frames.append(imageio.imread("temp_frame.png"))
    plt.clf()

# Salva i frame come animazione GIF
imageio.mimsave("animation.gif", frames, duration=0.05)
os.remove("temp_frame.png")
