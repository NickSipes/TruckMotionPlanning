class node:
	def __init__(self, x, y, cost, pind):
		self.x = x 			# x index
		self.y = y			# y index
		self.cost = cost	# cost
		self.pind = pind	# parent index

	def getX(self):
		return self.x

	def getY(self):
		return self.y

	def getCost(self):
		return self.cost

	def getParentIndex(self):
		return self.pind

	def setCost(self, new):
		self.cost = new

	def setParentIndex(self, new):
		self.pind = new