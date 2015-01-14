#
# A module for defining vectors in terms of simple lists. contain general vector methods
#

from math import *

class vector3:

	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z
		self.magnitude = self.mag()
	
	def __add__(self, rhs):
		return vector3(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)

	def __sub__(self, rhs):
		return vector3(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)

	def mag(self):
		return sqrt(pow(self.x, 2) + pow(self.y, 2) + pow(self.z, 2))
	
	def normalise(self):
		mag = self.mag()
		self.x = self.x/mag
		self.y = self.y/mag
		self.z = self.z/mag
	
	def scale(self, factor):
		self.x = self.x * factor
		self.y = self.y * factor
		self.z = self.z * factor
	
	def set_pos(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z

def vec3_dot_prod(vec31, vec32):
	out = vec31.x * vec32.x + vec31.y * vec32.y + vec31.z * vec32.z
	return out
 
def vec3_cross_prod(vec31, vec32):
	vec3out = vector3(0, 0, 0)
	vec3out.x = vec31.y * vec32.z - vec31.z * vec32.y
	vec3out.y = vec31.z * vec32.x - vec31.x * vec32.z
	vec3out.z = vec31.x * vec32.y - vec31.y * vec32.x
	return vec3out

def vec3_scale(vec3, scale_factor):
	vec3out = vector3(0, 0, 0)
	vec3out.x = vec3.x * scale_factor
	vec3out.y = vec3.y * scale_factor
	vec3out.z = vec3.z * scale_factor
	return vec3out
