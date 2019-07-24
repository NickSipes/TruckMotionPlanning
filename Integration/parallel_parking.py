
import trailer_hybrid_a_star
import matplotlib
import math
import time
matplotlib.use('GTK3Agg')
from matplotlib import pyplot as plt

# include("../src/trailer_hybrid_a_star.jl")

def main():
    # print ("start!!", PROGRAM_FILE)
    print("Start!!")

    # initial state
    sx = -15.0  # [m]
    sy = 8.0  # [m]
    syaw0 = math.radians(180.0)
    syaw1 = math.radians(180.0)

    # goal state
    gx = 1.0  # [m]
    gy = 3.0  # [m]
    gyaw0 = math.radians(0.0)
    gyaw1 = math.radians(0.0)

    # set obstacles
    ox = []
    oy = []

    for i in range(-25, 25):
        ox.append(i)
        oy.append(15.0)

    for i in range(-10, 10):
        ox.append(i)
        oy.append(0.0)

    for i in range(-25, -10):
        ox.append(i)
        oy.append(5.0)

    for i in range(10, 25):
        ox.append(i)
        oy.append(5.0)

    for i in range(0, 5):
        ox.append(10.0)
        oy.append(i)

    for i in range(0, 5):
        ox.append(-10.0)
        oy.append(i)


    oox = ox[:]
    ooy = oy[:]
    # plot(oox, ooy, ".k")
    # axis("equal")
    # show()

    # path generation
    # nor sure what to do with @time, so took it out from path below
    path = trailer_hybrid_a_star.calc_hybrid_astar_path(sx, sy, syaw0, syaw1, gx, gy, gyaw0, gyaw1, ox, oy, trailer_hybrid_a_star.XY_GRID_RESOLUTION, trailer_hybrid_a_star.YAW_GRID_RESOLUTION)

    # ====Animation=====
    show_animation(path, oox, ooy, sx, sy, syaw0, syaw1, gx, gy, gyaw0, gyaw1)

    print("Done!!")
    # print("Done!!", PROGRAM_FILE)

main()


def show_animation(path, oox, ooy, sx, sy, syaw0, syaw1, gx, gy, gyaw0, gyaw1):
    plt.plot(oox, ooy, ".k")
    trailer_hybrid_a_star.trailerlib.plot_trailer(sx, sy, syaw0, syaw1, 0.0)
    trailer_hybrid_a_star.trailerlib.plot_trailer(gx, gy, gyaw0, gyaw1, 0.0)
    x = path.x
    y = path.y
    yaw = path.yaw
    yaw1 = path.yaw1
    direction = path.direction

    steer = 0.0
    for ii in range(1, len(x)):
        plt.cla()
        plt.plot(oox, ooy, ".k")
        plt.plot(x, y, "-r", label="Hybrid A* path")

        if ii < (len(x) - 1):
            k = (yaw[ii + 1] - yaw[ii]) / trailer_hybrid_a_star.MOTION_RESOLUTION
            if not direction[ii]:
                k *= -1
            steer = math.atan2(trailer_hybrid_a_star.WB * k, 1.0)
        else:
            steer = 0.0


        trailer_hybrid_a_star.trailerlib.plot_trailer(x[ii], y[ii], yaw[ii], yaw1[ii], steer)
        plt.grid(True)
        plt.axis("equal")
        time.sleep(0.0001)

