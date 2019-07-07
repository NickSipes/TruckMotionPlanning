"""This file contains the majority of the functions for forming the path taken by the tractor trailer.  It is called
by a script that defines the start and goal state and outputs the path (with time) taken."""

import matplotlib as plt
import math
import numpy as np
import kdtree as kd

# TODO resolve the following Julia packages: DataFrames, NearestNeighbors (should be kdtree)

# TODO download each file from github to incorporate
import rs_path
import grid_a_star
import trailerlib

# Resolution and other variables
XY_GRID_RESOLUTION = 2                  # meters
YAW_GRID_RESOLUTION = math.radians(15)  # radians
GOAL_TYAW_TH = math.radians(5)          # radians
MOTION_RESOLUTION = 0.1                 # meters (path interpolation resolution)
N_STEER = 20                            # number of steer commands
EXTEND_AREA = 5                         # meters, used in calc_config to extend the configurations space beyond the
                                        # outermost obstacles
SKIP_COLLISION_CHECK = 4                # skip number for collision checks

# Cost variables
SB_COST = 100           # switchback cost (makes changing directions significantly expensive)
BACK_COST = 5           # backwards penalty cost (makes the algorithm favor going forward)
STEER_CHANGE_COST = 5   # cost to change steering angle
STEER_COST = 1          # cost to be turning
JACKKNIF_COST = 200     # cost if a jack knife occurs
H_COST = 5              # heuristic cost weighting

# Trailer specific parameters from trailerlib
WB = trailerlib.WB                # meters, wheel base
LT = trailerlib.LT                # meters, length of trailer
MAX_STEER = trailerlib.MAX_STEER  # radians, maximum steering angle


class Node: # TODO: implement which inputs can be left blank until needed
    def __init__(self, x_index, y_index, yaw_index, direction, x, y, yaw, yaw1, directions, steer, cost, parent_index):
        self.xind = x_index         # index of the x position (related to grid structure of Cspace) TODO: understand index
        self.yind = y_index         # index of the y position (related to grid structure of Cspace)
        self.yawind = yaw_index     # index of the yaw angle
        self.direction = direction  # direction of travel, true=forward, false=backward
        self.x = x                  # meters, x position
        self.y = y                  # meters, y position
        self.yaw = yaw              # radians, yaw angle of tractor
        self.yaw1 = yaw1            # radians, yaw angle of trailer
        self.directions = directions  # directions of positions and angles (TODO understand "directions")
        self.steer = steer          # steer input
        self.cost = cost            # cost TODO: which cost?
        self.pind = parent_index    # index of the parent node


class Config:
    """
    Class representing the Cspace bounds and resolution
    """
    def __init__(self, minx, miny, minyaw, minyawt, maxx, maxy, maxyaw, maxyawt, xw, yw, yaww, yawtw, xyreso, yawreso):
        # TODO understand each of these variables
        # TODO implement which inputs can be left blank until needed
        self.minx = minx
        self.miny = miny
        self.minyaw = minyaw
        self.minyawt = minyawt
        self.maxx = maxx
        self.maxy = maxy
        self.maxyaw = maxyaw
        self.maxyawt = maxyawt
        self.xw = xw
        self.yw = yw
        self.yaww = yaww
        self.yawtw = yawtw
        self.xyreso = xyreso
        self.yawreso = yawreso

class Path:
    def __init__(self, x, y, yaw, yaw1, direction, cost):
        self.x = x                  # meters, x position
        self.y = y                  # meters, y position
        self.yaw = yaw              # radians, tractor angle
        self.yaw1 = yaw1            # radians, trailer angle
        self.direction = direction  # direction of motion (true = forward, false = backward)
        self.cost = cost            # cost TODO which cost?

def calc_hybrid_astar_path(sx, sy, syaw, syaw1, gx, gy, gyaw, gyaw1, ox, oy, xyreso, yawreso):
    """
    TODO describe what the calc_hybrid_astar_path function does
    :param sx:      meters, start x position
    :param sy:      meters, start y position
    :param syaw:    radians, start tractor angle
    :param syaw1:   radians, start trailer angle
    :param gx:      meters, goal x position
    :param gy:      meters, goal y position
    :param gyaw:    radians, goal tractor angle
    :param gyaw1:   radians, goal trailer angle
    :param ox:      meters, list of x coordinates of each obstacle
    :param oy:      meters, list of y coordinates of each obstacle
    :param xyreso:  meters, xy grid resolution
    :param yawreso: radians, yaw resolution
    :return:
    """

    syaw, gyaw = rs_path.pi_2_pi(syaw), rs_path.pi_2_pi(gyaw) # function used to ensure yaw is between -pi and pi
    oxy = [(ox[i], oy[i]) for i in range(0, len(ox))]   # forms a list of tuples from the two lists ox and oy for
                                                        # input to kd.create

    # forms a kd tree from the xy points of the obstacles
    kdtree = kd.create(oxy)

    # Form the configuration space
    c = calc_config(ox, oy, xyreso, yawreso)

    # Define the starting node
    nstart = Node(round(sx / xyreso), round(sy / xyreso), round(syaw / yawreso), True, [sx], [sy], [syaw], [syaw1],
                  [True], 0, 0, -1)

    # Define the goal node
    ngoal = Node(round(gx / xyreso), round(gy / xyreso), round(gyaw / yawreso), True, [gx], [gy], [gyaw], [gyaw1],
                 [True], 0, 0, -1)

    # Determine holonomic costs of each xy coordinate on the grid (numpy dependent)
    h_dp = calc_holonomic_with_obstacle_heuristics(ngoal, ox, oy, xyreso)

    open = {}    # Dictionary to hold open nodes
    closed = {}  # Dictionary to hold closed nodes
    pq = {}  # Dictionary to hold queue of indexes and costs, used to determine which node in open to expand next
    fnode = []   # TODO determine if list is correct structure and update comment with purpose, Julia code was "nothing"

    # Add the start node to the dictionary of open nodes
    open[calc_index(nstart, c)] = nstart

    # Add the start node to the dictionary of indexes and costs
    pq[calc_index(nstart,c)] = calc_cost(nstart, h_dp, c)

    u, d = calc_motion_inputs()
    nmotion = len(u)


def calc_config(ox, oy, xyreso, yawreso):
    """
    This function returns an object of the config class that defines the bounds of the Cspace and the size of the grid
    on the Cspace based on the resolution of each dimension.
    :param ox: meters, list of x coordinates of each obstacle
    :param oy: meters, list of y coordinates of each obstacle
    :param xyreso: meters, xy grid resolution
    :param yawreso: radians, yaw resolution
    :return: config, class Config, the Cspace
    """
    # Expand configuration space beyond the outermost obstacles
    min_x_m = min(ox) - EXTEND_AREA
    min_y_m = min(oy) - EXTEND_AREA
    max_x_m = max(ox) + EXTEND_AREA
    max_y_m = max(oy) + EXTEND_AREA

    # TODO determine if ox and oy need to be updated globally
    # adds the two min and max xy points to the list of obstacles
    ox.append(min_x_m)
    oy.append(min_y_m)
    ox.append(max_x_m)
    oy.append(max_y_m)

    # used to form the size of the xy grid of the Cspace based on xy resolution
    minx = round(min_x_m/xyreso)
    miny = round(min_y_m/xyreso)
    maxx = round(max_x_m/xyreso)
    maxy = round(max_y_m/xyreso)

    # size of the xy grid of the Cspace
    xw = round(maxx - minx)
    yw = round(maxy - miny)

    # size of the yaw dimension of the Cspace based on yaw resolution
    minyaw = round(-math.pi/yawreso) - 1
    maxyaw = round(math.pi/yawreso)
    yaww = round(maxyaw - minyaw)

    # TODO determine what these t variants of yaw do
    minyawt = minyaw
    maxyawt = maxyaw
    yawtw = yaww

    # return the Cspace as an object of class config detailing the bounds of the Cspace and the resolution of each
    # dimension
    config = Config(minx, miny, minyaw, minyawt, maxx, maxy, maxyaw, maxyawt, xw, yw, yaww, yawtw, xyreso, yawreso)
    return config


def calc_holonomic_with_obstacle_heuristics(gnode, ox, oy, xyreso):
    """
    This function calculates the holonomic costs based on the function in grid_a_star
    :param gnode: Node, goal node
    :param ox: meters, list of x coordinates of each obstacle
    :param oy: meters, list of y coordinates of each obstacle
    :param xyreso: meters, xy grid resolution
    :return: 2D array with heuristic costs for each xy coordinate TODO ensure matches function output from grid_a_star
    """
    h_dp = grid_a_star.calc_dist_policy(gnode.x[-1], gnode.y[-1], ox, oy, xyreso, 1)  # TODO ensure inputs are correct format
    return h_dp


def calc_index(node, c):
    """
    This determines the key for the given node as a function of the index of each dimension
    :param node: Node, a node
    :param c: Config, the Cspace
    :return: integer, the key for the given node
    """
    # Account for the tractor position (x,y) and yaw
    ind = (node.yawind - c.minyaw)*c.xw*c.yw + (node.yind - c.miny)*c.xw + (node.xind - c.minx)

    # Account for the trailer yaw
    yaw1ind = round(node.yaw1[-1]/c.yawreso)
    ind += (yaw1ind - c.minyawt)*c.xw*c.yw*c.yaww

    if ind <= 0:
        print(f"Error: calculated index was invalid: {ind}")

    return ind


def calc_cost(n, h_dp, c):
    """
    This function calculates the sum of the past cost and heuristic cost
    :param n: Node, current node
    :param h_dp: 2D array, heuristic costs
    :param c: Config, Cspace
    :return: total cost of the current node
    """
    return n.cost + H_COST*h_dp[n.xind - c.minx, n.yind - c.miny]


def calc_motion_inputs():
    """
    This function calculates the possible inputs
    :return: u - list of steering inputs, d - list of drive inputs (only two speeds, forward/backward)
    """
    # Form list of all steering inputs
    steer = [i for i in np.arange(MAX_STEER/N_STEER, MAX_STEER, MAX_STEER/N_STEER)]

    # initialize steering input with straight
    u = [0]

    # initialize drive input with empty list
    d = []

    # Put all the +/- steer inputs into u
    for i in steer:
        u.append(i)
        u.append(-i)

    # for all the steering inputs form d for forward
    for i in range(len(u)):
        d.append(1)

    # for all the steering inputs form d for backward
    for i in range(len(u)):
        d.append(-1)

    # append u to u so when indexed, each steering input pairs with both a forward and back driving input
    for i in u:
        u.append(i)

    return u, d
