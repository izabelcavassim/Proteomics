from sys import argv 
import skimage
from skimage import io
import os
from scipy import ndimage #find the centroid
import numpy as np #matrix dealing
from skimage.draw import circle # to make the red spots
from skimage.color import label2rgb # to coloring
from skimage.morphology import disk #smoothing 
from skimage.filters import roberts # edge operators to find spots
from skimage.filters.rank import mean #smoothing the image
from skimage.filters.rank import gradient  #edge dectection method
import heapq # queue priority Q
from skimage import morphology # reducing complexity

skimage.data_dir = os.getcwd()
filename = os.path.join(skimage.data_dir, argv[1])
gel = io.imread(filename)

class PriorityQueue:
    def __init__(self):
        self._queue = []
        self._index = 0
        self._history = set()

    def push(self, item, priority):
        heapq.heappush(self._queue, (-priority, self._index, item))
        self._index += 1
        self._history.add(item)

    def pop(self):
        return heapq.heappop(self._queue)[-1]

# construct a matrix of 1's
def watershed(gel, method):
	matrix_one = np.ones((3,3)) 	
	M = np.zeros(gel.shape) #and a matrix of zeros 
	#making the gradient of the image
	if method == ' ':
		gel_gradient = skimage.filters.rank.gradient(gel, matrix_one)

	if method == 'roberts':
		gel_gradient = skimage.filters.roberts(gel)			

	#finding the coordinates of the dark spots:
	x, y = np.where(gel_gradient == gel_gradient.min()) #lowest points =darker
	coordinates = zip(x,y)
	# adding these points to queue with priority X[i,j]:
	q = PriorityQueue()
	counter = 0
	for i in coordinates:
		q.push(i, gel_gradient[i]) #inserting
		counter += 1
		M[i] = counter #labelling the matrix

# As long as the queue is not empty, search and store in N the elements around the point: (pop),
	while len(q._queue) > 0:
		last = q.pop()
		h, w = M.shape
		N = list()
		B =[(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1), (1,0), (1,1)] # relative coordinates 
		for bx, by in B:
			nx = last[0] + bx
			ny = last[1] + by
			if nx < 0 or nx >= h or ny < 0 or ny >= w: #introducing the exceptions
				continue
			if M[nx, ny] > 0:
				N.append(M[nx, ny])
			elif (nx,ny) not in q._history:
				coo = (nx,ny)
				q.push(coo, gel_gradient[coo])
	#see if the list N is unique
		if len(set(N)) == 1: 
			M[last[0], last[1]] = N[0]
	return M

M = watershed(gel, ' ')
print len(np.unique(M))

#improving the matrix M to reduce the noise:
print "Removed connected components smaller than 64 pixels."
R = morphology.remove_small_objects(M.astype(int), 64)
print 'The number of spots found for connectivity < 64 is:'
print len(np.unique(R))

unique = np.unique(R)
centroid = ndimage.measurements.center_of_mass(gel, R, unique)
# finding the intensity and priting results:
for i in range(len(centroid)):
	print 'Spot %d x: %f y: %f intensity: %d' % (i+1, centroid[i][0], centroid[i][1], gel[int(centroid[i][0]), int(centroid[i][1])])

colored = label2rgb(R, gel, bg_label=0)
for i in range(len(centroid)):
	rr, cc = circle(centroid[i][0], centroid[i][1], 2)
	colored[rr, cc] = (1,0,0)

io.imsave(argv[2], colored)

print "Smoothing the image before running watershed..."
loc_mean = mean(gel, disk(1))
smooth_M = watershed(loc_mean, ' ')
print "Number of spots found on the smoothed image without any post process is:"
print len(np.unique(smooth_M))
R = morphology.remove_small_objects(smooth_M.astype(int), 6)
colored = label2rgb(R, gel, bg_label=0)
io.imsave('smooth_cleaned.png', colored)
print "Number of spots found on the smoothed image after image processing is:"
print len(np.unique(R))

#############################
#Answering the question 1, using different edge operator methods to detect markers:
print 'The number of spots found using roberts gradient method is:'
M_roberts = watershed(gel, 'roberts')
R_roberts = morphology.remove_small_objects(M_roberts.astype(int), 64)
print len(np.unique(R_roberts))
colored = label2rgb(R_roberts, gel, bg_label=0)
io.imsave('roberts.png', colored)