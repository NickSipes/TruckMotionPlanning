"""This file contains the majority of the functions for forming the path taken by the tractor trailer.  It is called
by a script that defines the start and goal state and outputs the path (with time) taken."""

# import matplotlib as plt
import math
import numpy as np
import operator
from scipy import spatial


import rs_path
import Astar_Tractor_Trailer
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
JACKKNIF_COST = 200     # cost used to discourage jackknifes by making paths with hidh difference between yaw and yaw1
# undesirable
H_COST = 5              # heuristic cost weighting

# Trailer specific parameters from trailerlib
WB = trailerlib.WB                # meters, wheel base
LT = trailerlib.LT                # meters, length of trailer
MAX_STEER = trailerlib.MAX_STEER  # radians, maximum steering angle


class Node:
    def __init__(self, x_index, y_index, yaw_index, direction, x, y, yaw, yaw1, directions, steer, cost, parent_index):
        self.xind = x_index         # index of the x position (related to grid structure of Cspace)
        self.yind = y_index         # index of the y position (related to grid structure of Cspace)
        self.yawind = yaw_index     # index of the yaw angle
        self.direction = direction  # direction of travel, true=forward, false=backward
        self.x = x                  # meters, x position
        self.y = y                  # meters, y position
        self.yaw = yaw              # radians, yaw angle of tractor
        self.yaw1 = yaw1            # radians, yaw angle of trailer
        self.directions = directions  # list, direction associated with each x and y position
        self.steer = steer          # steer input
        self.cost = cost            # cost
        self.pind = parent_index    # index of the parent node


class Config:
    """
    Class representing the Cspace bounds and resolution
    """
    def __init__(self, minx, miny, minyaw, minyawt, maxx, maxy, maxyaw, maxyawt, xw, yw, yaww, yawtw, xyreso, yawreso):
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
        self.cost = cost            # cost


def calc_hybrid_astar_path(sx, sy, syaw, syaw1, gx, gy, gyaw, gyaw1, ox, oy, xyreso, yawreso):
    """

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
    kdtree = spatial.KDTree(oxy)

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

    opened = {}    # Dictionary to hold open nodes
    closed = {}  # Dictionary to hold closed nodes
    pq = {}  # Dictionary to hold queue of indexes and costs, used to determine which node in open to expand next
    fnode = []

    # Add the start node to the dictionary of open nodes
    opened[calc_index(nstart, c)] = nstart

    # Add the start node to the dictionary of indexes and costs
    pq[calc_index(nstart,c)] = calc_cost(nstart, h_dp, c)

    u, d = calc_motion_inputs()
    nmotion = len(u)

    while True:

        # Exit while loop and return nothing if open is expended and a solution is not found
        if not opened:
            print("Error: Cannot find path, No open nodes remaining.")
            return []

        # Obtain the index of the next node
        c_id = min(pq.keys())  # TODO fixed syntax

        # Removed the obtained index from the queue
        pq.pop(c_id)

        # Get the node with the obtained index and remove it from open
        current = opened.get(c_id)
        opened.pop(c_id)

        # Add the current node to the closed list
        closed.update({c_id: current})

        # get full data of current node and isupdated flag
        isupdated, fpath = update_node_with_analystic_expantion(current, ngoal, c, ox, oy, kdtree, gyaw1)

        if isupdated: # goal has been found
            fnode = fpath
            break  # exit the while loop

        inityaw1 = current.yaw1[0]  # TODO fixed index error

        # cycle through all inputs and check if each next node is in closed, else add it to open if not already there
        for i in range(nmotion):
            node = calc_next_node(current, c_id, u[i], d[i], c)
            if not verify_index(node, c, ox, oy, inityaw1, kdtree):
                continue  # continue with next node since this node's configuration is invalid

            node_ind = calc_index(node, c)

            # Check if node is already in the closed dictionary
            if node_ind in closed:
                continue  # continue with next node since this one is already closed

            # Check if node is already in the open dictionary
            if node_ind not in opened:
                opened.update({node_ind: node})  # add the node to open
                pq.update({node_ind: calc_cost(node, h_dp, c)})  # add the node cost information to pq
            else:
                if opened[node_ind].cost > node.cost:  # if node is in the open dictionary, but the new node cost is
                    # lower, update the node in open
                    opened[node_ind] = node

    print(f"number of nodes in open and closed = {len(opened) + len(closed)}")

    path = get_final_path(closed, fnode, nstart, c)

    return path


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
    :return: 2D array with heuristic costs for each xy coordinate
    """
    h_dp = Astar_Tractor_Trailer.calculate_dist_policy(gnode.x[-1], gnode.y[-1], ox, oy, xyreso, 1)
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
    return n.cost + H_COST*h_dp[n.xind - c.minx][n.yind - c.miny]  # TODO fixed indexing syntax


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
    u.extend(u)  # TODO fixed (old syntax was an infinite loop)

    return u, d


def update_node_with_analystic_expantion(current, ngoal, c, ox, oy, kdtree, gyaw1):
    """
    This function updates the current node with "analystic" data
    :param current: Node, current node
    :param ngoal: Node, goal node
    :param c: Config, configuration space
    :param ox: meters, list of x positions of each obstacle
    :param oy: meters, list of y positions of each obstacle
    :param kdtree: kdtree, KD tree made from ox and oy
    :param gyaw1: radians, goal position trailer yaw
    :return:
    """

    apath = analystic_expantion(current, ngoal, c, ox, oy, kdtree)

    # Formulate data for the "f" node
    if apath:  # if apath = [], skip the if statement since there is no path (ignore the unresolved attribute reference
        # errors as they won't happen since apath is only a list if it is empty and then the if statement is false and
        # the block doesn't execute)
        fx = apath.x[2:-1]
        fy = apath.y[2:-1]
        fyaw = apath.yaw[2:-1]
        steps = [MOTION_RESOLUTION*direct for direct in apath.directions]  # TODO fixed element wise multiplication
        yaw1 = trailerlib.calc_trailer_yaw_from_xyyaw(apath.x, apath.y, apath.yaw, current.yaw1[-1], steps)
        # check if the trailer yaw is outside the trailer yaw threshold
        if abs(rs_path.pi_2_pi(yaw1[-1] - gyaw1)) >= GOAL_TYAW_TH:
            return False, []  # current node is not the goal node based on trailer alone
        fcost = current.cost + calc_rs_path_cost(apath, yaw1)  # cost to next node
        fyaw1 = yaw1[2:-1] # array of trailer yaws
        fpind = calc_index(current, c)  # index of parent node
        fd = []  # array of directions
        for d in apath.directions[2:-1]:
            if d >= 0:
                fd.append(True)
            else:
                fd.append(False)

        fsteer = 0

        fpath = Node(current.xind, current.yind, current.yawind, current.direction, fx, fy, fyaw, fyaw1, fd, fsteer,
                     fcost, fpind)

        return True, fpath

    return False, []


def analystic_expantion(n, ngoal, c, ox, oy, kdtree):
    """
    This function determines the least costly path to the next node that doesn't collide with an obstacle
    :param n: Node, current node
    :param ngoal: Node, goal node
    :param c: Config, Cspace
    :param ox: meters, list of x positions of obstacles
    :param oy: meters, list of y positions of obstacles
    :param kdtree: kdtree, KD Tree
    :return: an RS path leading to the next least costly node that doesn't collide with an obstacle
    """

    # Extract the x,y, and yaw info from the current node
    sx = n.x[-1]
    sy = n.y[-1]
    syaw = n.yaw[-1]

    # Calculate the sharpest turn possible
    max_curvature = math.tan(MAX_STEER)/WB

    # Call rs_path to determine the paths available based on current and goal node
    paths = rs_path.calc_paths(sx, sy, syaw, ngoal.x[-1], ngoal.y[-1], ngoal.yaw[-1], max_curvature,
                               step_size=MOTION_RESOLUTION)

    # If there is no available path, stop and return nothing
    if not paths:
        return []

    pathqueue = {}  # Dictionary for holding the path and yaw1 information

    # Put the rs path and yaw1 combinations into pathqueue
    for path in paths:
        steps = [MOTION_RESOLUTION*direct for direct in path.directions]  # TODO fixed element wise multiplication
        yaw1 = trailerlib.calc_trailer_yaw_from_xyyaw(path.x, path.y, path.yaw, n.yaw1[-1], steps)
        pathqueue.update({path: calc_rs_path_cost(path, yaw1)})

    # Go through each path, starting with the lowest cost, and check for collisions, return the first viable path
    for i in range(len(pathqueue)):
        for path_key, path_value in pathqueue.items():  # extract path with lowest cost TODO changed since old syntax did not work for a dictionary (this line + 3 are new)
            if path_value == min(pathqueue.values()):
                path = path_key
                break
        pathqueue.pop(path)  # remove the lowest cost path from the queue

        steps = [MOTION_RESOLUTION*direct for direct in path.directions]  # TODO fixed element wise multiplication
        yaw1 = trailerlib.calc_trailer_yaw_from_xyyaw(path.x, path.y, path.yaw, n.yaw1[-1], steps)
        ind = [i for i in np.arange(1, len(path.x), SKIP_COLLISION_CHECK)]
        path_x = [path.x[i] for i in ind]  # TODO added this and following 4 lines to correctly form list inputs to the collision checker
        path_y = [path.y[i] for i in ind]
        path_yaw = [path.yaw[i] for i in ind]
        path_yaw1 = [yaw1[i] for i in ind]
        if trailerlib.check_trailer_collision(ox, oy, path_x, path_y, path_yaw, path_yaw1, kdtree=kdtree):
            return path

    return []


def calc_rs_path_cost(rspath, yaw1):
    """
    This function calculates the total cost for a given RS path and yaw1
    :param rspath: rs_path.Path, RS path
    :param yaw1: radians, trailer yaw
    :return: total cost
    """

    # Initialize cost
    cost = 0

    # Update cost for length and direction of path
    for l in rspath.lengths:
        if l >= 0:
            cost += l
        else:
            cost += abs(l)*BACK_COST

    # Update cost if direction changes in the path (detected by change in sign of adjacent lengths)
    for i in (range(len(rspath.lengths) - 1)):
        if rspath.lengths[i]*rspath.lengths[i+1] < 0:
            cost += SB_COST

    # Update cost for curved path
    for ctype in rspath.ctypes:
        if ctype != "S":
            cost += STEER_COST*abs(MAX_STEER)

    # Form list of steering inputs
    nctypes = len(rspath.ctypes)  # number of ctypes in the given path
    ulist = []  # list used to convert the string ctypes to numbers for cost calculation
    for i in range(nctypes):
        if rspath.ctypes[i] == "R":
            ulist.append(-MAX_STEER)
        elif rspath.ctypes[i] == "L":
            ulist.append(MAX_STEER)
        else:  # TODO I added this else statement because ulist was initialized as empty vice an array of zeros
            ulist.append(0)

    # Update cost for changing direction of turn
    for i in range(nctypes - 2):  # TODO fixed max index value (changed from -1 to -2)
        cost += STEER_CHANGE_COST*abs(ulist[i+1] - ulist[i])

    # Update cost to prevent jackknifes (the most the trailer folds toward the truck, the higher the cost)
    cost += JACKKNIF_COST*sum(np.absolute(rs_path.pi_2_pi(rspath.yaw-yaw1)))  # TODO changed to np.absolute to allow element wise operation

    return cost


def calc_next_node(current, c_id, u, d, c):
    """
    This function calculates the next node based on the current node and the given inputs.
    :param current: Node, current node
    :param c_id: index, index of current node
    :param u: input, steering input
    :param d: input, driving input
    :param c: Config, Cspace
    :return:
    """

    arc_l = XY_GRID_RESOLUTION*1.5

    nlist = math.floor(arc_l/MOTION_RESOLUTION) + 1

    # Initial motion from the current node
    xlist = [current.x[-1] + d*MOTION_RESOLUTION*math.cos(current.yaw[-1])]
    ylist = [current.y[-1] + d*MOTION_RESOLUTION*math.sin(current.yaw[-1])]
    yawlist = [rs_path.pi_2_pi(current.yaw[-1] + d*MOTION_RESOLUTION*math.tan(u)/WB)]
    yaw1list = [rs_path.pi_2_pi(current.yaw1[-1] + d*MOTION_RESOLUTION*math.sin(current.yaw[-1] - current.yaw1[-1])/LT)]

    # Discrete path from current node to the next node for the given inputs
    for i in range(nlist-1):
        xlist.append(xlist[-1] + d*MOTION_RESOLUTION*math.cos(yawlist[-1]))
        ylist.append(ylist[-1] + d*MOTION_RESOLUTION*math.sin(yawlist[-1]))
        yawlist.append(rs_path.pi_2_pi(yawlist[-1] + d*MOTION_RESOLUTION*math.tan(u)/WB))
        yaw1list.append(rs_path.pi_2_pi(yaw1list[-1] + d*MOTION_RESOLUTION*math.sin(yawlist[-1] - yaw1list[-1])/LT))

    # Determine the index in each dimension of the next node
    xind = round(xlist[-1]/c.xyreso)
    yind = round(ylist[-1]/c.xyreso)
    yawind = round(yawlist[-1]/c.yawreso)

    # Calculate cost to the next node
    addedcost = 0  # Initialize cost of moving from current node to the next node

    # Cost of traveling length in given direction
    if d > 0:
        direction = True
        addedcost += abs(arc_l)
    else:
        direction = False
        addedcost += abs(arc_l)*BACK_COST

    # Cost of switching direction
    if direction != current.direction:  # A direction change occured
        addedcost += SB_COST

    # Cost of steering
    addedcost += STEER_COST*abs(u)

    # Cost of changing direction
    addedcost += STEER_CHANGE_COST*abs(current.steer - u)

    # Cost of acute angle between trailer and tractor
    addedcost += JACKKNIF_COST*sum(np.absolute(rs_path.pi_2_pi(np.subtract(yawlist, yaw1list))))  # TODO changed to np.absolute and np.subtract to allow element wise operation

    cost = current.cost + addedcost

    # Form list of directions for the node class
    directions = [direction for i in range(len(xlist))]

    node = Node(xind, yind, yawind, direction, xlist, ylist, yawlist, yaw1list, directions, u, cost, c_id)

    return node


def verify_index(node, c, ox, oy, inityaw1, kdtree):
    """
    This function verifies that the x and y index of the given node are valid and collision free
    :param node: Node, node to check x and y index of
    :param c: Config, Cspace
    :param ox: meters, list of x positions of obstacles
    :param oy: meters, list of y positions of obstacles
    :param inityaw1: radians, first increment of yaw1 from the current node configuration
    :param kdtree: scipy.spatial.KDTree, KD Tree of obstacles
    :return: boolean, validity of the x and y indexes
    """

    # Check if the x index is valid
    if (node.xind - c.minx) >= c.xw:  # x index is too high
        return False
    elif (node.xind - c.minx) <= 0:   # x index is too low
        return False

    if (node.yind - c.miny) >= c.yw:  # y index is too high
        return False
    elif (node.yind - c.miny) <= 0:   # y index is too low
        return False

    # Check if the node collides with an obstacle
    steps = [MOTION_RESOLUTION*direct for direct in node.directions]  # TODO fixed element wise multiplication
    yaw1 = trailerlib.calc_trailer_yaw_from_xyyaw(node.x, node.y, node.yaw, inityaw1, steps)
    ind = [i for i in np.arange(1, len(node.x), SKIP_COLLISION_CHECK)]
    node_x = [node.x[i] for i in ind]  # TODO added this and following 4 lines to correctly form list inputs to the collision checker
    node_y = [node.y[i] for i in ind]
    node_yaw = [node.yaw[i] for i in ind]
    node_yaw1 = [yaw1[i] for i in ind]
    if not trailerlib.check_trailer_collision(ox, oy, node_x, node_y, node_yaw, node_yaw1, kdtree=kdtree):
        return False

    # If none of the above returned false, return true
    return True


def get_final_path(closed, ngoal, nstart, c):
    """
    Form the path from the start to the goal based on the closed list (possible due to parent index in each node)
    :param closed: dict of Node, all closed nodes
    :param ngoal: Node, goal node
    :param nstart: Node, start node
    :param c: Config, Cspace
    :return:
    """

    # Initialize recursive list of path parameters
    rx = ngoal.x.copy()  # TODO edited reverse statements to correctly save the array to rx, ry, ryaw, and directions
    rx.reverse()
    ry = ngoal.y.copy()
    ry.reverse()
    ryaw = ngoal.yaw.copy()
    ryaw.reverse()
    ryaw1 = np.flip(ngoal.yaw1)  # TODO fixed "reverse" since yaw1 is a nparray, not list
    direction = ngoal.directions.copy()
    direction.reverse()
    nid = ngoal.pind
    finalcost = ngoal.cost

    # Form recursive list of path parameters
    while True:
        n = closed[nid]
        rx_temp = n.x.copy()  # TODO Fixed reverse and appened statements for rx, ry, ryaw, ryaw1, and direction
        rx_temp.reverse()
        rx.extend(rx_temp)
        ry_temp = n.y.copy()
        ry_temp.reverse()
        ry.extend(ry_temp)
        ryaw_temp = n.yaw.copy()
        ryaw_temp.reverse()
        ryaw.extend(ryaw_temp)
        np.append(ryaw1, np.flip(n.yaw1))
        direction_temp = n.directions.copy()
        direction_temp.reverse()
        direction.extend(direction_temp)
        nid = n.pind
        if is_same_grid(n, nstart):
            break

    # Reverse the path parameters to start at start and end at goal
    rx.reverse()
    ry.reverse()
    ryaw.reverse()
    np.flip(ryaw1)
    direction.reverse()

    direction[0] = direction[1]

    path = Path(rx, ry, ryaw, ryaw1, direction, finalcost)

    return path


def is_same_grid(node1, node2):
    """
    This function determines if two given nodes are in the same grid based on the index (based on resolution)
    :param node1: Node, a node
    :param node2: Node, another node
    :return: boolean, status of being in the same grid
    """

    if node1.xind != node2.xind:
        return False

    if node1.yind != node2.yind:
        return False

    if node1.yawind != node2.yawind:
        return False

    return True
