from sys import argv 
import skimage
import os 
from skimage.measure import label #labeling 
from scipy import ndimage #to find the centroid 
from skimage import io
from skimage.color import label2rgb # coloring
import numpy as np # to deal with matrix
from skimage.draw import circle # to make the red spots 
import matplotlib.pyplot as plt #ploting

skimage.data_dir = os.getcwd()
filename = os.path.join(skimage.data_dir, argv[1])
gel = io.imread(filename)

#defining the threshold:
thresh = gel.mean() - gel.std()*3 # 3 is k 

#a matrix of true or false
label1 = label(gel <= thresh)

#coloring the image by given the threshold points
colored = label2rgb(label1, gel, bg_label=0) # use these labels and give them a color

#calculating the number of spots
num_spot = np.unique(label1).max()
print 'The number of spots detected is:'
print num_spot

centroid = ndimage.measurements.center_of_mass(gel, label1, range(1, num_spot +1))

# finding the intensity and priting results:
for i in range(len(centroid)):
	print 'Spot %d x: %f y: %f intensity: %d' % (i+1, centroid[i][0], centroid[i][1], gel[int(centroid[i][0]), int(centroid[i][1])])

#coloring the the center of the spots with red color
for i in range(len(centroid)):
	rr, cc = circle(centroid[i][0], centroid[i][1], 2)
	colored[rr, cc] = (1,0,0)

io.imsave(argv[2], colored)

print '##############################'
############## answering the question 
# question 1
def threshold(gel, thresh):
	label1 = label(gel <= thresh)
	num_spot = np.unique(label1).max()
	return num_spot

threshs = [0,50,100,150,200]
spots = []
for i in threshs:
	spots.append(threshold(gel, i))
print 'The number of spots found with the changing of the threshold'
print spots

#question 2
plt.plot(threshs, spots, 'ro')
plt.ylabel('Spots')
plt.xlabel('Threshold')
plt.grid(True)
plt.show()

#question 3
ks = [1,2,3,4,5]
thre = []
labels = []
num_spots = []
for k in ks:
	thre.append(gel.mean() - gel.std()*k)
	labels.append(label(gel <= thre[k-1]))
	num_spots.append(np.unique(labels[k-1]).max())
print 'The number of spots found with the changing of k' 
print num_spots

#question 4
#adding noise in the matrix:
# Draw samples from a Poisson distribution.
def adding_noise(lam, matrix):
	noise = np.random.poisson(lam, matrix.shape).astype(float)
	image = matrix + noise
	num_spot = np.unique(image).max()
	return np.unique(image).max()
noises = [1, 3, 5, 10, 30, 50]
for n in noises:
	print 'Number of spots found for the given lambda %.1f' % (n)
	print adding_noise(n, label1)

#As lambda increases, the Poisson distribution approaches a normal distribution. And make the 
# image more clear.