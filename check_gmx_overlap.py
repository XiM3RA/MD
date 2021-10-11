#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import glob

plt.ion()
class Graph():
	legendlist = []
	def __call__(self):
		"""
		Invoked in the main loop, opens initial display window
		and reads .dat files into list
		"""
		self.directories = glob.glob("*dat")
		plt.show()

	def pf(self, sel):
		"""
		Updates the graph by reading file of choice
		"""
		self.x_axis = []
		self.y_axis = []
		self.legendlist.append(self.directories[sel].rstrip(".dat"))
		with open (self.directories[sel], 'r') as fh:
			for line in fh:
				temp = line.split()
				self.x_axis.append(float(temp[0]))
				self.y_axis.append(float(temp[1]))
		plt.plot(self.x_axis, self.y_axis)
		plt.legend(self.legendlist, loc='upper right')

	def ls(self):
		"""
		Analagous to UNIX ls command
		"""
		for z, item in enumerate(self.directories):
			print(z, item)
	def clear(self):
		plt.clf()

g = Graph()
g()
a = "Commands for usage\ng.ls() -- list all the data files available\ng.pf({number}) -- give this function a number, get a plot\ng.clear() -- clears display"
print(a)
