'''cs410 camera calibration assignment
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import axes3d
from numpy.linalg import eig

def calibrateCamera3d(data):
	"""Calculates perspective projection matrix for given data"""
	# Create an ampty numpy matrix to store the measurement matrix, A
	N = data.shape[0]
	A = np.empty([(2*(N)),12])
	# Iterate through the data matrix and set the values of the measurement matrix
	for i in range(N):
		index = 2*i 					# Each row of data = two rows of A
		P = np.matrix(data[i,0:3])		# A1 -> Xi Yi Zi
		P = np.hstack([P, [[1]]])		# A1 -> Xi Yi Zi 1
		zeros = np.zeros((1,4))			
		lin = np.hstack([P, zeros])		# A1 -> Xi Yi Zi 1 0 0 0 0
		# obtain xi and yi to perform multiplications
		x = data[i,3]				
		y = data[i,4]
		xEnd = np.multiply(-x,P)		# -xiXi -xiYi -xiZi -xi
		yEnd = np.multiply(-y,P)		# -yiXi -yiYi -yiZi -yi
		lin = np.hstack([lin, xEnd])	
		A[index] = lin 					# A1 -> Xi Yi Zi 1 0 0 0 0 -xiXi -xiYi -xiZi -xi
		index = index + 1 				# next row of A using same data row
		lin2 = np.hstack([zeros,P,yEnd])
		A[index] = lin2					# A2 -> 0 0 0 0 Xi Yi Zi 1 -yiXi -yiYi -yiZi -yi
	AtA =  ((A.transpose()).dot(A))		# Transpose of A dot A for 12 x 12 matrix
	d,v = eig(AtA)						# Eigenvalues and eigenvectors of AtA
	v = v[:,-1]/v[1,-1]					# Last eigenvalue is minimum, take corresponding vector
	M = np.vstack([v[0:4],v[4:8],v[8:12]])
	# Stack the first, middle and last four elements for perspective projection matrix, M
	return M
		

def visualiseCameraCalibration3D(data, P):
	"""Renders a 2D plot showing 
	i) the measured 2D image point and 
	ii) the reprojection of the 3D calibraion points as computed by P """
	coords = matrixP(data)	# Matrix of world coordinates
	p = P.dot(coords) 		# p = mP from notes
	w = p[-1]				# The last element of each row, w
	p = p / w[None,:]		# Dived each column by w to scale
	fig = plt.figure()
	ax = fig.gca()
	ax.set_title('Measured Image Points and Reprojected Points Computed by P')
	ax.axis([0, 800, 0, 700])
	ax.plot(data[:,3], data[:,4],'r.')
	ax.plot(p[0,:], p[1,:], 'b.')
	blue_dots = mpatches.Patch(color='blue', label='Reprojected points')
	red_dots = mpatches.Patch(color='red', label='Measured points')
	plt.legend(handles=[red_dots, blue_dots])
	plt.show()

def evalutateCameraCalibration3D(data, P):
	"""Prints the mean, variance, minimum and maximum distances in 
	pixels between the measured and reprojected image feature locations"""
	coords = matrixP(data)  # Matrix of world coordinates
	p = P.dot(coords)		# p = mP from notes
	w = p[-1]				# The last element of each row, w
	p = p / w[None,:]		# Dived each column by w to scale

	# Create empty numpy matrix to store the error of each point
	errors = np.empty([p.shape[1],1])
	for i in range(p.shape[1]):
	# The error of each point is the hypothenuse of the distance in x and y
		x = p[0,i] - data[i,3]
		y = p[1,i] - data[i,4]
		errors[i] = (x**2 + y**2)**0.5

	print ' Mean error \t= ', np.mean(errors), '\tpixels'
	print ' Variance \t= ', np.var(errors), '\tpixels'
	print ' Minimum error \t= ', np.amin(errors), '\tpixels'
	print ' Maximum error \t= ', np.amax(errors), '\tpixels'

def matrixP(data):
	"""Creates 'P' matrix of the World Coords given in data"""
	P = data[:,0:3] 				# X Y Z of each point
	P = np.insert(P,3,1, axis=1)	# -> X Y Z 1 
	P = P.transpose()				# Each of above bceoms a column (required for later multiplication)
	return P

if __name__ == "__main__":
	print "\nStarting programme...\n"
	data = np.loadtxt('data.txt')
	M = calibrateCamera3d(data)
	evalutateCameraCalibration3D(data,M)
	visualiseCameraCalibration3D(data, M)