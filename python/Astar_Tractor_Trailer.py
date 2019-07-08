#
# Grid based A* shortest path planning
#
# RBE 550: Motion Planning
# 

from matplotlib import pyplot
import math
from scipy import spatial
from nodeObject import node

def calculate_dist_policy(goalX, goalY, obstacleX, obstacleY, gridResolution, vehicleRadius):

	goalNode = node(goalX / gridResolution, goalY / gridResolution, 0, -1)

	# Using a list
	obstacleX = [x/gridResolution for x in obstacleX]
	obstacleY = [y/gridResolution for y in obstacleY]

	obstacleMap, minX, minY, maxX, maxY, xWidth, yWidth = calculate_obstacle_map(obstacleX, obstacleY, gridResolution, vr)

	# Create dictionaries of open and closed nodes
	openList = dict()
	closedList = dict()

	# Put the goal node on the openList
	openList.insert({calculate_node_index(goalNode, xWidth, minX, minY) : goalNode})

	# Get the motion : cost association of our model
	motion = get_motion_model()
	numberOfMotions = length(motion)

	# Prioirity Queue (Not Sure What it does)
	pq = dict({calculate_node_index(goalNode, xWidth, minX, minY): goalNode.getCost()})

	# Meat of the method
	while True:
		if len(openList) == 0:
			print('Search is done')
			break

		# Set current node to one with maximum key value
		currentKey 	= max(pq.keys())
		currentNode = openList.pop(currentKey)

		# Expand search grid based on motion model
		for i in range(numberOfMotions):

			# Create the next node after executing the motion in the ith row
			searchNode = node(currentNode.getX() + motion[i,1], currentNode.getY() + motion[i,2], currentNode.getCost() + motion[i,3], currentKey)

			# If the node is not verified, check if the next motion will produce a valid motion
			if not verify_node(searchNode, minX, minY, xWidth, yWidth, obstacleMap):
				continue

			currentNodeKey = calculate_node_index(searchNode, xWidth, yWidth, minX, minY)

			# If the node has already been explored, skip over it to the next motion
			if currentNodeKey in closedList.keys():
				continue

			# If the node is in open, see if cost can be updated
			if currentNodeKey in openList.keys():

				# If a cheaper path has been found through the open node
				if openList[currentNodeKey].getCost() > searchNode.getCost():

					# Update with the new parent
					openList[currentNodeKey].setCost(searchNode.getCost())
					openList[currentNodeKey].setParentIndex(currentKey)
		
			else:
				# Add to open set
				openList[currentNodeKey] = node

				# Add to queue with 
				pq.update({calculate_node_index(searchNode, xWidth, minX, minY) : searchNode.getCost()})

	# Calculate the policy map using the closed nodes
	pmap = calc_policy_map(closedList, xWidth, yWidth, minX, minY)

	return pmap


def calc_policy_map(closedList, xWidth, yWidth, minX, minY):

	# Create an array of dimension xWidth by Ywidth filled with values of infinity
	pmap = [[math.inf for j in range(xWidth)] for i in range(yWidth)]

	# Iterate through the dictionary items in the closed list
	for value in closedList.values():
		pmap[value.getX() - minX, value.getY() - minY] = value.getCost()

	return pmap

def calc_astar_path(startX, startY, goalX, goalY, obstacleX, obstacleY, gridResolution, vr):

	# Set all of the input values to the approriate resolution of our grid
	# *****
	startNode = node(startX/gridResolution, startY/gridResolution, 0, -1)
	goalNode = node(goalX/gridResolution, goalY/gridResolution, 0, -1)

	obstacleX = [x/gridResolution for x in obstacleX]
	obstacleY = [y/gridResolution for y in obstacleY]

	# Get the obstacle map object with the dimensions
	obstacleMap, minX, minY, maxX, maxY, xWidth, yWidth = calc_obstacle_map(obstacleX, obstacleY, gridResolution, vr)

	# Open and closed set
	openList = dict()
	closedList = dict()

	# Insert the start node into the open list
	openList.update({calc_index(startNode, xWidth, minX, minY) : startNode})

	motion = get_motion_model()
	numberOfMotions = len(motion)

	# Queue the start node with a cost determined by the heuristic distance between start and goal node
	pq = dict({calc_index(startNode, xWidth, minX, minY): calc_cost(startNode, goalNode)})

	# Meat of the method
	while True:
		if len(openList) == 0:
			print('ERROR: No open set')
			break

		# Set current node to one with maximum key value and remove from open list
		currentKey 	= max(pq.keys())
		pq.pop(currentKey)
		currentNode = openList.pop(currentKey)

		# Check if we are at the goal node
		if currentNode.getX() == goalNode.getX() and currentNode.getY() == goalNode.getY():
			print("GOAL!!")

			# Insert the goal node into the closed list
			closedList.update({currentKey : currentNode})
			break

		# Add to closed list
		closedList.update({currentKey : currentNode})

		# Expand search grid based on motion model
		for i in range(numberOfMotions):

			# Create the next node after executing the motion in the ith row
			searchNode = node(currentNode.getX() + motion[i][0], currentNode.getY() + motion[i][1], currentNode.getCost() + motion[i][2], currentKey)

			# If the node is not verified, check if the next motion will produce a valid motion
			if not verify_node(searchNode, minX, minY, xWidth, yWidth, obstacleMap):
				# The "continue" statement terminates the current iteration of the loop and procees to the next
				continue

			# Calculate the index for our current node
			searchNodeKey = calc_index(searchNode, xWidth, minX, minY)

			# If the node has already been explored, skip over it to the next motion
			if searchNodeKey in closedList.keys():
				continue

			# If the node is in open, see if cost can be updated
			if searchNodeKey in openList.keys():

				# If a cheaper path has been found through the open node
				if openList[searchNodeKey].getCost() > searchNode.getCost():
					# Update with the new parent
					openList[searchNodeKey].setCost(searchNode.getCost())
					openList[searchNodeKey].setParentIndex(currentKey)
		
			else:
				# Add to open set
				openList.update({searchNodeKey : searchNode})

				# Add to queue with 
				pq.update({searchNodeKey : calc_cost(searchNode, goalNode)})

	# Calculate the policy map using the closed nodes
	pathX, pathY = get_final_path(closedList, goalNode, startNode, xWidth, minX, minY, gridResolution)

	return pathX, pathY


# Check that the node is a navigable space
def verify_node(node, minX, minY, xWidth, yWidth, obstacleMap):

	# Check if the node is outside map bounds
	# ***
	if node.getX() - minX >= xWidth:
		return False
	elif node.getX() - minX <= 0:
		return False

	if node.getY() - minY >= yWidth:
		return False
	elif node.getY() - minY <= 0:
		return False
	# ***

	# Check if the node is colliding with an obstacle
	if obstacleMap[int(node.getX() - minX)][int(node.getY() - minY)]:
		return False

	# Node has been successfuly verified
	return True

# Calculate the cost of the node
def calc_cost(node, goalNode):
	return (node.getCost() + heuristic_cost(node.getX() - goalNode.getX(), node.getY() - goalNode.getY()))


# Return the motion model of the tractor-trailer
def get_motion_model():

	# Return a 3 X 8 matrix that associates a unique movement with a particular cost
	# dx, dy, cost
	motion = [[1, 0, 1],
			[0, 1, 1],
			[-1, 0, 1],
			[0, -1, 1],
			[-1, -1, math.sqrt(2)],
			[-1, 1, math.sqrt(2)],
			[1, -1, math.sqrt(2)],
			[1, 1, math.sqrt(2)]]

	return motion

# Calculate the index of a node based on map size
def calc_index(node, xWidth, minX, minY):
	return (node.getY() - minY) * xWidth + (node.getX() - minX)

def calc_obstacle_map(obstacleX, obstacleY, gridResolution, vehicleRadius):

	# Get obstacle bounds
	minX = min(obstacleX)
	minY = min(obstacleY)
	maxX = max(obstacleX)
	maxY = max(obstacleY)

	xWidth = maxX - minX
	yWidth = maxY - minY

	# Create a 2D list where each element is False
	obstacleMap = [[False for j in range(int(xWidth) + 1)] for i in range(int(yWidth) + 1)]

	# Create a KDTree with the x and y positions of the obstacles on the map
	# zip() relates values of similar index position in both lists to be a pair
	# list() converts the result into a list for KDTree()
	kdtree = spatial.KDTree(list(zip(obstacleX, obstacleY)))

	for ix in range(int(xWidth) + 1):
		x = ix + minX

		for iy in range(int(yWidth) + 1):
			y = iy + minY

			# Return the position of the nearest neighbor to x,y and the distance between them
			distance, indexPosition =  kdtree.query([x, y])

			# If the nearest neighbor is closer than the size of the vehicle, accounting for grid resolution,
			# make it an obstacle
			if distance <= vehicleRadius / gridResolution:
				obstacleMap[ix][iy] = True

	return obstacleMap, minX, minY, maxX, maxY, xWidth, yWidth

def get_final_path(closedList, goalNode, startNode, xWidth, minX, minY, gridResolution):

	pathX = [goalNode.getX()]
	pathY = [goalNode.getY()]
	nodeID = calc_index(goalNode, xWidth, minX, minY)

	while True:

		# Find the node on the closed list with the associated ID
		n = closedList[nodeID]

		# Add the x and y positions of the node to the final path
		pathX.append(n.getX())
		pathY.append(n.getY())

		# Update the ID with the current nodes predecessor
		nodeID = n.getParentIndex()

		# If the last node added is the start node, terminate
		if pathX[len(pathX) - 1] == startNode.getX() and pathY[len(pathY) - 1] == startNode.getY():
			print('PATH FOUND')
			break

	# Set the nodes to the scale of the grid resoultion
	pathX = [gridResolution * n for n in pathX]
	pathY = [gridResolution * n for n in pathY]

	# Reverse the order so the first point is the start Node
	pathX.reverse()
	pathY.reverse()

	return pathX, pathY

# Find the node, within the open queue, with the lowest cost
def search_min_cost_node(open, goalNode):

	# Initalize the tracking variables
	minNode = 0
	minCost = math.inf

	# Search all nodes in the open list
	for n in open:
		temp_cost = n.getCost() + heuristic(n.getX() - goalNode.getX(), n.getY() - goalNode.getY())

		# If n node is cheaper than any node so far, save it
		if minCost > cost:
			minNode = n
			minCost = cost

	return minNode	

# Calculate heuristic cost for the node
def heuristic_cost(x,y):
	return math.sqrt(math.pow(x , 2) + math.pow(y, 2))

def main():

	startX = 10.0 # meters
	startY = 10.0 # meters
	goalX  = 50.0 # meters
	goalY  = 50.0 # meters 

	obstacleX = []
	obstacleY = []

	# Build the map of obstacles.
	# Each x,y pair is an obstacle
	for i in range(60):
		obstacleX.append(i)
		obstacleY.append(0.0)

	for i in range(60):
		obstacleX.append(60.0)
		obstacleY.append(i)

	for i in range(60):
		obstacleX.append(i)
		obstacleY.append(60.0)

	for i in range(60):
		obstacleX.append(0.0)
		obstacleY.append(i)

	for i in range(40):
		obstacleX.append(20.0)
		obstacleY.append(i)

	for i in range(40):
		obstacleX.append(40.0)
		obstacleY.append(60.0 - i)

	VEHICLE_RADIUS  = 5.0 # meters
	GRID_RESOLUTION = 1.0 # meters

	# Compute the path using the start, goal, obstacle location, grid resolution, and vehicle radius
	pathX, pathY = calc_astar_path(startX, startY, goalX, goalY, obstacleX, obstacleY, GRID_RESOLUTION, VEHICLE_RADIUS)

	fig = pyplot.figure()

	# Plot obstacles, start, goal, and A* path
	pyplot.scatter(obstacleX, obstacleY)
	pyplot.scatter(pathX, pathY)
	pyplot.plot( startX, startY, 'g', goalX, goalY, 'c')

	# Show the final plots
	pyplot.show()


# THIS IS WHAT RUNS WHEN THE FILE IS RAN
main()




















