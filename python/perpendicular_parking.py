# This function describes the environment for a perpendicular parking scenario and executes it.

import trailer_hybrid_a_star
import matplotlib
import math
import time
matplotlib.use('GTK3Agg')
from matplotlib import pyplot as plt


def show_animation(path, oox, ooy, sx, sy, syaw0, syaw1, gx, gy, gyaw0, gyaw1):
    """
    Function to animate the tractor trailer along the path found by the algorithm.
    :param path: Path object, output from the trailer_hybrid_a_star.calc_hybrid_astar_path function
    :param oox: list, x coordinates of obstacles
    :param ooy: list, y coordinates of obstacles
    :param sx: float, starting x coordinate
    :param sy: float, starting y coordinate
    :param syaw0: float, starting tractor angle
    :param syaw1: float, starting trailer angle
    :param gx: float, goal x coordinate
    :param gy: float, goal y coordinate
    :param gyaw0: float, goal tractor angle
    :param gyaw1: float, goal traier angle
    :return: a pretty animation
    """
    # Plot the obstacles, start, and goal
    plt.plot(oox, ooy, ".k")
    trailer_hybrid_a_star.trailerlib.plot_trailer(sx, sy, syaw0, syaw1, 0)
    trailer_hybrid_a_star.trailerlib.plot_trailer(gx, gy, gyaw0, gyaw1, 0)

    # Form the lists containing the path information
    x = path.x
    y = path.y
    yaw = path.yaw
    yaw1 = path.yaw1
    direction = path.direction

    steer = 0

    # Plot the path
    for ii in range(x):
        plt.cla()
        plt.plot(oox, ooy, ".k")
        plt.plot(x, y, "-r", label="Hybrid A* path")

        if ii < len(x)-1:  # Not the last point in the path
            k = (yaw[ii + 1] - yaw[ii])/trailer_hybrid_a_star.MOTION_RESOLUTION
            if not direction[ii]:  # If going backward
                k*= -1
            steer = math.atan2(trailer_hybrid_a_star.WB*k, 1)
        else:
            steer = 0
        trailer_hybrid_a_star.trailerlib.plot_trailer(x[ii], y[ii], yaw[ii], yaw1[ii], steer)
        plt.grid(True)
        plt.axis("equal")
        time.sleep(0.0001)


def main():
    """
    Function to fully execute a given scenario.
    :return:
    """
    # Initial State
    sx = 14  # Starting x (m)
    sy = 10  # Starting y (m)
    syaw0 = math.radians(0)  # Starting tractor angle (rad)
    syaw1 = math.radians(0)  # Starting trailer angle (rad)

    # Goal State
    gx = 0  # Goal x (m)
    gy = 0  # Goal y (m)
    gyaw0 = math.radians(0)  # Goal tractor angle (rad)
    gyaw1 = math.radians(0)  # Goal trailer angle (rad)

    # Obstacles

    ox = []  # x coordinates of obstacles
    oy = []  # y coordinates of obstacles

    for i in range(-25, 25):  # Wall from (-25, 15) to (25, 15)
        ox.append(i)
        oy.append(15)

    for i in range(-25, -4):  # Wall from (-25, 4) to (-4, 4)
        ox.append(i)
        oy.append(4)

    for i in range(-15, 4):  # Wall from (-4, -15) to (-4, 4)
        ox.append(-4)
        oy.append(i)

    for i in range(-15, 4):  # Wall from (4, -15) to (4, 4)
        ox.append(4)
        oy.append(i)

    for i in range(4, 25):  # Wall from (4, 4) to (4, 25)
        ox.append(i)
        oy.append(4)

    for i in range(-4, 4):  # Wall from (-4, -15) to (4, -15)
        ox.append(i)
        oy.append(15)

    # oox = ox[:]  # TODO are these lines necessary, they seem to just copy the ox and oy lists
    # ooy = oy[:]

    # Generate Path
    # TODO @time is a timing function.  That can be added later with the time module.
    path = trailer_hybrid_a_star.calc_hybrid_astar_path(sx, sy, syaw0, syaw1, gx, gy, gyaw0, gyaw1, ox, oy, trailer_hybrid_a_star.XY_GRID_RESOLUTION, trailer_hybrid_a_star.YAW_GRID_RESOLUTION)

    # Animate Path
    show_animation(path, ox, oy, sx, sy, syaw0, syaw1, gx, gy, gyaw0, gyaw1)

    # Executed code
    main()
