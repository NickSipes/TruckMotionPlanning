import math
import numpy as np
import matplotlib
matplotlib.use('GTK3Agg')
from matplotlib import pyplot as plt

STEP_SIZE = 0.1
MAX_PATH_LENGTH = 1000

class Path:
    def __init__(self, lengths, ctypes, L, x, y, yaw, directions):
        self.lengths = lengths
        self.ctypes = ctypes
        self.L = L
        self.x = x
        self.y = y
        self.yaw = yaw
        self. directions = directions

# Utility Functions
def set_path(paths, lengths, ctypes):
    """
    This function adds a path with the given lengths and curves if not already in the list of paths
    :param paths: list of Path, list of paths
    :param lengths: list of lengths forming a path
    :param ctypes: list of curves forming a path
    :return: list of Paths, paths
    """
    path = Path(lengths, ctypes, 0, [], [], [], [])

    # Check if path already in paths
    for tpath in paths:
        type_is_same = (tpath.ctypes == path.ctypes)
        if type_is_same:
            if sum(tlength - length for (tlength, length) in zip(tpath.lengths, path.lengths)) <= 0.01:
                return paths  # no need to add path

    path.L = sum([abs(i) for i in lengths])

    if path.L >= MAX_PATH_LENGTH:
        return paths

    if path.L < 0.01:
        print(f"Error: Path length is {path.L}")

    paths.append(path)

    return paths


def mod2pi(phi):
    """
    Function to limit the given angle to +/- pi
    :param phi: radians, angle to limit
    :return: radians, limited angle
    """
    v = phi % 2*math.pi
    if v < -math.pi:
        v += 2*math.pi
    elif v > math.pi:
        v -= 2*math.pi
    return v


def cart2pol(x, y):
    r = math.sqrt(x**2 + y**2)
    t = math.atan2(y, x)
    return r, t


def calc_tau_omega(u, v, xi, eta, phi):
    """

    :param u:
    :param v:
    :param xi:
    :param eta:
    :param phi:
    :return:
    """
    delta = mod2pi((u - v))
    a = math.sin(u) - math.sin(delta)
    b = math.cos(u) - math.cos(delta)- 1

    t1 = math.atan2(eta*a - xi*b, xi*a + eta*b)
    t2 = 2*(math.cos(delta) - math.cos(v) - math.cos(u)) + 3

    if t2 < 0:
        tau = mod2pi(t1 + math.pi)
    else:
        tau = mod2pi(t1)

    omega = mod2pi(tau - u + v - phi)

    return tau, omega


#  TODO changed entire pi_2_pi function to handle lists
def pi_2_pi(angles):
    """

    :param angles: float or list of floats
    :return:
    """
    if isinstance(angles, float):
        while angles > math.pi:
            angles -= 2*math.pi

        while angles < -math.pi:
            angles += 2*math.pi
        return angles
    else:  # input was a list
        angle_out = []
        for iangle in angles:
            while iangle > math.pi:
                iangle -= 2*math.pi

            while iangle < -math.pi:
                iangle += 2*math.pi
            angle_out.append(iangle)
        return angle_out


def get_label(path):
    """

    :param path:
    :return:
    """
    label = ""

    for (m, l) in zip(path.ctypes, path.lengths):
        label = label + m
        if l > 0:
            label = label + "+"
        else:
            label = label + "-"

    return label


# Path Shape Functions
def SLS(x, y, phi):
    """

    :param x:
    :param y:
    :param phi:
    :return:
    """
    phi = mod2pi(phi)
    if y > 0 and phi > 0 and phi < math.pi*0.99:
        xd = x - y/math.tan(phi)
        t = xd - math.tan(phi/2)
        u = phi
        v = math.sqrt((x - xd)**2 + y**2) - math.tan(phi/2)
        return True, t, u, v
    elif y < 0 and phi > 0 and phi < math.pi*0.99:
        xd = x - y / math.tan(phi)
        t = xd - math.tan(phi / 2)
        u = phi
        v = -math.sqrt((x - xd) ** 2 + y ** 2) - math.tan(phi / 2)
        return True, t, u, v
    return False, 0, 0, 0


def LSL(x, y, phi):
    """

    :param x:
    :param y:
    :param phi:
    :return:
    """
    u, t = cart2pol(x - math.sin(phi), y - 1 + math.cos(phi))
    if t >= 0:
        v = mod2pi(phi - t)
        if v >= 0:
            return True, t, u, v

    return False, 0, 0, 0


def LRL(x, y, phi):
    """

    :param x:
    :param y:
    :param phi:
    :return:
    """
    u1, t1 = cart2pol(x - math.sin(phi), y - 1 + math.cos(phi))

    if u1 <= 4:
        u = -2*math.asin(0.25*u1)
        t = mod2pi(t1 + 0.5*u + math.pi)
        v = mod2pi(phi - t + u)

        if t >= 0 and u <= 0:
            return True, t, u, v

    return False, 0, 0, 0


def LRLRn(x, y, phi):
    """

    :param x:
    :param y:
    :param phi:
    :return:
    """
    xi = x + math.sin(phi)
    eta = y - 1 - math.cos(phi)
    rho = 0.25*(2 + math.sqrt(xi**2 + eta**2))

    if rho <= 1:
        u = math.acos(rho)
        t, v = calc_tau_omega(u, -u, xi, eta, phi)
        if t >= 0 and v <= 0:
            return True, t, u, v

    return False, 0, 0, 0


def LRLRp(x, y, phi):
    """

    :param x:
    :param y:
    :param phi:
    :return:
    """
    xi = x + math.sin(phi)
    eta = y - 1 - math.cos(phi)
    rho = (20 - xi**2 - eta**2)/16

    if 0 <= rho <= 1:
        u = -math.acos(rho)
        if u >= -0.5*math.pi:
            t, v = calc_tau_omega(u, u, xi, eta, phi)
            if t >= 0 and v <= 0:
                return True, t, u, v

    return False, 0, 0, 0


def LRSL(x, y, phi):
    """

    :param x:
    :param y:
    :param phi:
    :return:
    """
    xi = x - math.sin(phi)
    eta = y - 1 + math.cos(phi)
    rho, theta = cart2pol(xi, eta)

    if rho >= 2:
        r = math.sqrt(rho**2 - 4)
        u = 2 - r
        t = mod2pi(theta + math.atan2(r, -2))
        v = mod2pi(phi - 0.5*math.pi - t)
        if t >= 0 and u <= 0 and v >= 0:
            return True, t, u, v

    return False, 0, 0, 0


def LRSR(x, y, phi):
    """

    :param x:
    :param y:
    :param phi:
    :return:
    """
    xi = x + math.sin(phi)
    eta = y - 1 - math.cos(phi)
    rho, theta = cart2pol(-eta, xi)

    if rho >= 2:
        t = theta
        u = 2 - rho
        v = mod2pi(t + 0.5*math.pi - phi)
        if t >= 0 and u <= 0 and v <= 0:
            return True, t, u, v

    return False, 0, 0, 0


def LRSLR(x, y, phi):
    """

    :param x:
    :param y:
    :param phi:
    :return:
    """
    xi = x + math.sin(phi)
    eta = y - 1 - math.cos(phi)
    rho, theta = cart2pol(xi, eta)
    if rho >= 2:
        u = 4 - math.sqrt(rho**2 - 4)
        if u <= 0:
            t = mod2pi(math.atan2((4-u)*xi - 2*eta, -2*xi + (u - 4)*eta))  # TODO fixed parenthesis error
            v = mod2pi(t - phi)

            if t >= 0 and v >= 0:
                return True, t, u, v

    return False, 0, 0, 0


def interpolate(ind, l, m, maxc, ox, oy, oyaw, px, py, pyaw, directions):
    """

    :param ind:
    :param l:
    :param m:
    :param maxc:
    :param ox:
    :param oy:
    :param oyaw:
    :param px:
    :param py:
    :param pyaw:
    :param directions:
    :return:
    """
    if m == "S":
        px[ind] = ox + l*math.cos(oyaw)/maxc
        py[ind] = oy + l*math.sin(oyaw)/maxc
        pyaw[ind] = oyaw
    else:
        ldx = math.sin(l)/maxc
        if m =="L":
            ldy = (1 - math.cos(l))/maxc
            pyaw[ind] = oyaw + l
        elif m == "R":
            ldy = (1 - math.cos(l))/(-maxc)
            pyaw[ind] = oyaw - l
        gdx = math.cos(-oyaw)*ldx + math.sin(-oyaw)*ldy
        gdy = -math.sin(-oyaw)*ldx + math.cos(-oyaw)*ldy
        px[ind] = ox + gdx
        py[ind] = oy + gdy

    if l > 0:
        directions[ind] = 1
    else:
        directions[ind] = -1

    return px, py, pyaw, directions


# Path Functions
def SCS(x, y, phi, paths):
    """

    :param x:
    :param y:
    :param phi:
    :param paths:
    :return:
    """
    flag, t, u, v = SLS(x, y, phi)

    if flag:
        paths = set_path(paths, [t, u, v], ["S", "L", "S"])  # TODO fixed lower case s

    flag, t, u, v = SLS(x, -y, -phi)

    if flag:
        paths = set_path(paths, [t, u, v], ["S", "R", "S"])

    return paths


def CSC(x, y, phi, paths):
    """

    :param x:
    :param y:
    :param phi:
    :param paths:
    :return:
    """
    flag, t, u, v = LSL(x, y, phi)
    if flag:
        paths = set_path(paths, [t, u, v], ["L", "S", "L"])

    flag, t, u, v = LSL(-x, y, -phi)
    if flag:
        paths = set_path(paths, [-t, -u, -v], ["L", "S", "L"])

    flag, t, u, v = LSL(x, -y, -phi)
    if flag:
        paths = set_path(paths, [t, u, v], ["R", "S", "R"])

    flag, t, u, v = LSL(-x, -y, phi)
    if flag:
        paths = set_path(paths, [-t, -u, -v], ["R", "S", "R"])

    flag, t, u, v = LSL(x, y, phi)
    if flag:
        paths = set_path(paths, [t, u, v], ["L", "S", "R"])

    flag, t, u, v = LSL(-x, y, -phi)
    if flag:
        paths = set_path(paths, [-t, -u, -v], ["L", "S", "R"])

    flag, t, u, v = LSL(x, -y, -phi)
    if flag:
        paths = set_path(paths, [t, u, v], ["R", "S", "L"])

    flag, t, u, v = LSL(-x, -y, phi)
    if flag:
        paths = set_path(paths, [-t, -u, -v], ["R", "S", "L"])

    return paths


def CCC(x, y, phi, paths):
    """

    :param x:
    :param y:
    :param phi:
    :param paths:
    :return:
    """
    flag, t, u, v = LRL(x, y, phi)
    if flag:
        paths = set_path(paths, [t, u, v], ["L", "R", "L"])

    flag, t, u, v = LRL(-x, y, -phi)
    if flag:
        paths = set_path(paths, [-t, -u, -v], ["L", "R", "L"])

    flag, t, u, v = LRL(x, -y, -phi)
    if flag:
        paths = set_path(paths, [t, u, v], ["R", "L", "R"])

    flag, t, u, v = LRL(-x, -y, phi)
    if flag:
        paths = set_path(paths, [-t, -u, -v], ["R", "L", "R"])

    # Backward direction
    xb = x*math.cos(phi) + y*math.sin(phi)
    yb = x*math.sin(phi) - y*math.cos(phi)

    flag, t, u, v = LRL(xb, yb, phi)
    if flag:
        paths = set_path(paths, [t, u, v], ["L", "R", "L"])

    flag, t, u, v = LRL(-xb, yb, -phi)
    if flag:
        paths = set_path(paths, [-t, -u, -v], ["L", "R", "L"])

    flag, t, u, v = LRL(xb, -yb, -phi)
    if flag:
        paths = set_path(paths, [t, u, v], ["R", "L", "R"])

    flag, t, u, v = LRL(-xb, -yb, phi)
    if flag:
        paths = set_path(paths, [-t, -u, -v], ["R", "L", "R"])

    return paths


def CCCC(x, y, phi, paths):
    """

    :param x:
    :param y:
    :param phi:
    :param paths:
    :return:
    """
    flag, t, u, v = LRLRn(x, y, phi)
    if flag:
        paths = set_path(paths, [t, u, -u, v], ["L", "R", "L", "R"])

    flag, t, u, v = LRLRn(-x, y, -phi)
    if flag:
        paths = set_path(paths, [-t, -u, u, -v], ["L", "R", "L", "R"])

    flag, t, u, v = LRLRn(x, -y, -phi)
    if flag:
        paths = set_path(paths, [t, u, -u, v], ["R", "L", "R", "L"])

    flag, t, u, v = LRLRn(-x, -y, phi)
    if flag:
        paths = set_path(paths, [-t, -u, u, -v], ["R", "L", "R", "L"])

    flag, t, u, v = LRLRp(x, y, phi)
    if flag:
        paths = set_path(paths, [t, u, -u, v], ["L", "R", "L", "R"])

    flag, t, u, v = LRLRp(-x, y, -phi)
    if flag:
        paths = set_path(paths, [-t, -u, u, -v], ["L", "R", "L", "R"])

    flag, t, u, v = LRLRp(x, -y, -phi)
    if flag:
        paths = set_path(paths, [t, u, -u, v], ["R", "L", "R", "L"])

    flag, t, u, v = LRLRp(-x, -y, phi)
    if flag:
        paths = set_path(paths, [-t, -u, u, -v], ["R", "L", "R", "L"])

    return paths


def CCSC(x, y, phi, paths):
    """

    :param x:
    :param y:
    :param phi:
    :param paths:
    :return:
    """
    flag, t, u, v = LRSL(x, y, phi)
    if flag:
        paths = set_path(paths, [t, -0.5 * math.pi, u, v], ["L", "R", "S", "L"])

    flag, t, u, v = LRSL(-x, y, -phi)
    if flag:
        paths = set_path(paths, [-t, 0.5 * math.pi, -u, -v], ["L", "R", "S", "L"])

    flag, t, u, v = LRSL(x, -y, -phi)
    if flag:
        paths = set_path(paths, [t, -0.5 * math.pi, u, v], ["R", "L", "S", "R"])

    flag, t, u, v = LRSL(-x, -y, phi)
    if flag:
        paths = set_path(paths, [-t, 0.5 * math.pi, -u, -v], ["R", "L", "S", "R"])

    flag, t, u, v = LRSR(x, y, phi)
    if flag:
        paths = set_path(paths, [t, -0.5 * math.pi, u, v], ["R", "L", "S", "R"])

    flag, t, u, v = LRSR(-x, y, -phi)
    if flag:
        paths = set_path(paths, [-t, 0.5 * math.pi, -u, -v], ["R", "L", "S", "R"])

    flag, t, u, v = LRSR(x, -y, -phi)
    if flag:
        paths = set_path(paths, [t, -0.5 * math.pi, u, v], ["L", "R", "S", "L"])

    flag, t, u, v = LRSR(-x, -y, phi)
    if flag:
        paths = set_path(paths, [-t, 0.5 * math.pi, -u, -v], ["L", "R", "S", "L"])

    # Backwards direction
    xb = x*math.cos(phi) + y*math.sin(phi)
    yb = x*math.sin(phi) - y*math.cos(phi)

    flag, t, u, v = LRSL(xb, yb, phi)
    if flag:
        paths = set_path(paths, [t, -0.5 * math.pi, u, v], ["L", "R", "S", "L"])

    flag, t, u, v = LRSL(-xb, yb, -phi)
    if flag:
        paths = set_path(paths, [-t, 0.5 * math.pi, -u, -v], ["L", "R", "S", "L"])

    flag, t, u, v = LRSL(xb, -yb, -phi)
    if flag:
        paths = set_path(paths, [t, -0.5 * math.pi, u, v], ["R", "L", "S", "R"])

    flag, t, u, v = LRSL(-xb, -yb, phi)
    if flag:
        paths = set_path(paths, [-t, 0.5 * math.pi, -u, -v], ["R", "L", "S", "R"])

    flag, t, u, v = LRSR(xb, yb, phi)
    if flag:
        paths = set_path(paths, [t, -0.5 * math.pi, u, v], ["R", "L", "S", "R"])

    flag, t, u, v = LRSR(-xb, yb, -phi)
    if flag:
        paths = set_path(paths, [-t, 0.5 * math.pi, -u, -v], ["R", "L", "S", "R"])

    flag, t, u, v = LRSR(xb, -yb, -phi)
    if flag:
        paths = set_path(paths, [t, -0.5 * math.pi, u, v], ["L", "R", "S", "L"])

    flag, t, u, v = LRSR(-xb, -yb, phi)
    if flag:
        paths = set_path(paths, [-t, 0.5 * math.pi, -u, -v], ["L", "R", "S", "L"])

    return paths

def CCSCC(x, y, phi, paths):
    """

    :param x:
    :param y:
    :param phi:
    :param paths:
    :return:
    """
    flag, t, u, v = LRSLR(x, y, phi)
    if flag:
        paths = set_path(paths, [t, -0.5*math.pi, u, -0.5*math.pi, v], ["L", "R", "S", "L", "R"])

    flag, t, u, v = LRSLR(-x, y, -phi)
    if flag:
        paths = set_path(paths, [-t, 0.5*math.pi, -u, 0.5*math.pi, -v], ["L", "R", "S", "L", "R"])

    flag, t, u, v = LRSLR(x, -y, -phi)
    if flag:
        paths = set_path(paths, [t, -0.5*math.pi, u, -0.5*math.pi, v], ["R", "L", "S", "R", "L"])

    flag, t, u, v = LRSLR(-x, -y, phi)
    if flag:
        paths = set_path(paths, [-t, 0.5*math.pi, -u, 0.5*math.pi, -v], ["R", "L", "S", "R", "L"])

    return paths


# High Level Path Functions
def generate_path(q0, q1, maxc):
    """

    :param q0:
    :param q1:
    :param maxc:
    :return:
    """
    # Calcuate delta in each dimension
    dx = q1[0] - q0[0]
    dy = q1[1] - q0[1]
    dth = q1[2] - q0[2]

    # Precalculate sin and cos
    c = math.cos(q0[2])
    s = math.sin(q0[2])
    x = (c*dx + s*dy)*maxc
    y = (-s*dx + c*dy)*maxc

    paths = []
    paths = SCS(x, y, dth, paths)
    paths = CSC(x, y, dth, paths)
    paths = CCC(x, y, dth, paths)
    paths = CCCC(x, y, dth, paths)
    paths = CCSC(x, y, dth, paths)
    paths = CCSCC(x, y, dth, paths)

    return paths


def generate_local_course(L, lengths, mode, maxc, step_size):
    """

    :param L:
    :param lengths:
    :param mode:
    :param maxc:
    :param step_size:
    :return:
    """
    npoint = math.floor(L/step_size) + len(lengths) + 3

    px = [0]*npoint
    py = [0]*npoint
    pyaw = [0]*npoint
    directions = [0]*npoint
    ind = 1  # Julia counts from 1, python counts from 0

    if lengths[0] > 0:
        directions[0] = 1
        d = step_size
    else:
        directions[0] = -1
        d = -step_size

    pd = d
    ll = 0

    for (m, l, i) in zip(mode, lengths, range(len(mode))):
        if l > 0:
            d = step_size
        else:
            d = -step_size

        # Set origin state
        ox, oy, oyaw = px[ind], py[ind], pyaw[ind]

        ind -= 1
        if i >= 2 and (lengths[i-1]*lengths[i]) > 0:
            pd = -d - ll
        else:
            pd = d - ll

        while abs(pd) <= abs(l):
            ind += 1
            px, py, pyaw, directions = interpolate(ind, pd, m, maxc, ox, oy, oyaw, px, py, pyaw, directions)
            pd += d

        ll = l - pd - d  # remaining length

        ind += 1

        px, py, pyaw, directions = interpolate(ind, l, m, maxc, ox, oy, oyaw, px, py, pyaw, directions)

    while px[-1] == 0:
        px.pop(-1)
        py.pop(-1)
        pyaw.pop(-1)
        directions.pop(-1)

    return px, py, pyaw, directions


def calc_paths(sx, sy, syaw, gx, gy, gyaw, maxc, step_size=STEP_SIZE):
    """

    :param sx:
    :param sy:
    :param syaw:
    :param gx:
    :param gy:
    :param gyaw:
    :param maxc:
    :param step_size:
    :return:
    """
    q0 = [sx, sy, syaw]
    q1 = [gx, gy, gyaw]

    paths = generate_path(q0, q1, maxc)
    for path in paths:
        x, y, yaw, directions = generate_local_course(path.L, path.lengths, path.ctypes, maxc, step_size*maxc)
        path.x = [math.cos(-q0[2])*ix + math.sin(-q0[2])*iy + q0[0] for (ix, iy) in zip(x, y)]
        path.y = [-math.sin(-q0[2])*ix + math.cos(-q0[2])*iy + q0[1] for (ix, iy) in zip(x, y)]
        path.yaw = pi_2_pi([iyaw + q0[2] for iyaw in yaw])
        path.directions = directions
        path.lengths = [l/maxc for l in path.lengths]
        path.L = path.L/maxc

    return paths


def calc_shortest_path(sx, sy,syaw, gx, gy, gyaw, maxc, step_size=STEP_SIZE):
    """

    :param sx:
    :param sy:
    :param syaw:
    :param gx:
    :param gy:
    :param gyaw:
    :param maxc:
    :param step_size:
    :return:
    """
    paths = calc_paths(sx, sy, syaw, gx, gy, gyaw, maxc, step_size=step_size)

    minL = math.inf
    best_path_index = -1
    for i in range(len(paths)):
        if paths[i].L <= minL:
            minL = paths[i].L
            best_path_index = i

    return paths[best_path_index]


def calc_curvature(x, y, yaw, directions):
    """

    :param x:
    :param y:
    :param yaw:
    :param directions:
    :return:
    """
    c = []
    ds = []

    for i in range(1,len(x)-1):
        dxn = x[i] - x[i-1]
        dxp = x[i+1] - x[i]
        dyn = y[i] - y[i-1]
        dyp = y[i+1] - y[i]
        dn = math.sqrt(dxn**2 + dyn**2)
        dp = math.sqrt(dxp**2 + dyp**2)
        dx = (dp*dxn/dn + dn*dxp/dp)/(dn + dp)
        ddx = 2*(dxp/dp - dxn/dn)/(dn + dp)
        dy = (dp*dyn/dn + dn*dyp/dp)/(dn + dp)
        ddy = 2*(dyp/dp - dyn/dn)/(dn + dp)
        curvature = (ddy*dx - ddx*dy)/(dx**2 + dy**2)
        d = (dn + dp)/2

        if np.isnan(curvature):
            curvature = 0

        if directions[i] <= 0:
            curvature = -curvature

        if len(c) == 0:
            c.append(curvature)
            ds.append(d)

        c.append(curvature)
        ds.append(d)

    ds.append(ds[-1])
    c.append(c[-1])

    return c, ds

"""

def main():
    # Initial and Final state
    start_x = 3
    start_y = 10
    start_yaw = math.radians(40)
    end_x = 0
    end_y = 1
    end_yaw = math.radians(0)
    max_curvature = 0.1

    # TODO @time just tracks the time required to execute and memory allocated (do we need that functionality?)
    bpath = calc_shortest_path(start_x, start_y, start_yaw, end_x, end_y, end_yaw, max_curvature)

    rc, rds = calc_curvature(bpath.x, bpath.y, bpath.yaw, bpath.directions)

    plt.subplot(121)
    plt.plot(bpath.x, bpath.y, "-r", label=get_label(bpath))

    plt.plot(start_x, start_y)
    plt.plot(end_x, end_y)

    plt.legend()
    plt.grid(True)
    plt.axis("equal")

    plt.subplot(122)
    plt.plot(rc, ".r", label="Reeds Shepp")
    plt.grid(True)
    plt.title("Curvature")

    plt.show()


main()
"""
