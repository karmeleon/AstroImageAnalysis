# Based on ASTR406 Informatics problem 2-3

# call as python analyze.py [fits image path] [standard deviations to search]

# this code is very slow, but it's readable and it works

import astropy.io.fits as fits
import numpy
import math
from scipy.interpolate import UnivariateSpline
from matplotlib import pyplot as pyplt
from matplotlib.colors import LogNorm as lm
from matplotlib import cm as cm

# finds the center of the star, as well as all the pixels part of the star
def center(image, sigmas):
	# find the star by finding pixels a given number of StD's above the mean brightness
	mean = numpy.mean(image)
	std = numpy.std(image)

	threshold = mean + std * sigmas

	x_pixels = []
	y_pixels = []

	ylen, xlen = image.shape

	xmean = ymean = totval = 0

	for i in range(xlen):
		for j in range(ylen):
			if image[j, i] > threshold:
				x_pixels.append(i)
				y_pixels.append(j)

				xmean += i * image[j, i]
				ymean += j * image[j, i]
				totval += image[j, i]


	minX = min(x_pixels)
	maxX = max(x_pixels)
	minY = min(y_pixels)
	maxY = max(y_pixels)

	width = max(maxX - minX, maxY - minY)

	# returns (xc, yc, x, y, window width, total)
	return (xmean / float(totval), ymean / float(totval), x_pixels, y_pixels, width, totval)

# returns flux inside radius about given center
def aperture(image, xc, yc, r):
	flux = 0
	ylen, xlen = image.shape

	for i in range(xlen):
		for j in range(ylen):
			if math.sqrt((i - xc)**2 + (j - yc)**2) <= r:
				flux += image[j, i]

	return flux

# returns pixels at 5, 25, 50, 75, and 90% of peak intensity
def isophotes(image, sigmas):
	peak = max([max(row) for row in image])

	threshold = numpy.mean(image) + numpy.std(image) * sigmas

	iso = []
	iso.append([(x, y) for (y, x), data in numpy.ndenumerate(image) if .05 * peak <= data < .25 * peak and data > threshold])
	iso.append([(x, y) for (y, x), data in numpy.ndenumerate(image) if .25 * peak <= data < .5 * peak and data > threshold])
	iso.append([(x, y) for (y, x), data in numpy.ndenumerate(image) if .5 * peak <= data < .75 * peak and data > threshold])
	iso.append([(x, y) for (y, x), data in numpy.ndenumerate(image) if .75 * peak <= data < .9 * peak and data > threshold])
	iso.append([(x, y) for (y, x), data in numpy.ndenumerate(image) if .9 * peak <= data and data > threshold])

	return iso

# returns a slice of the image along the x direction and the FWHM of the slice
def imageSlice(image, xc, yc, width):
	ylen, xlen = image.shape

	xstart = max(0, xc - width / 2)
	xend = min(xlen, xc + width / 2)

	# I cannot get array slicing to work for the life of me
	slice = []
	for i in range(int(xstart), int(xend)):
		slice.append(image[yc, i])

	# https://stackoverflow.com/questions/10582795/finding-the-full-width-half-maximum-of-a-peak/10583774#10583774

	shiftedSlice = []

	halfMax = max(slice) / 2
	baseline = numpy.mean(image)

	for y in slice:
		shiftedSlice.append(y - halfMax - baseline)

	x = numpy.linspace(0, width, width)
	spline = UnivariateSpline(x, shiftedSlice, s=0)
	r1, r2 = spline.roots()
	#r1 = 0
	#r2 = 0

	return (slice, r2 - r1)

def plot(image, xc, yc, width, slice, fwhm, iso, fluxes, radii, filename):
	file = ''
	if '/' in filename:
		file = filename[filename.rfind('/') + 1:]
	else:
		file = filename

	fig = pyplt.figure()
	fig.suptitle(file)
	sliceGraph = fig.add_subplot(2, 2, 1)
	imagePlot = fig.add_subplot(2, 2, 2)
	isophotePlot = fig.add_subplot(2, 2, 3)
	apertureCDF = fig.add_subplot(2, 2, 4)

	# slice graph
	sliceGraph.plot(numpy.linspace(0, len(slice), len(slice)), slice)
	sliceGraph.set_xlabel("Pixels")
	sliceGraph.set_ylabel("Count")

	sliceGraph.set_xlim(0, len(slice))
	sliceGraph.set_ylim(0, max(slice))

	sliceGraph.text(.05 * len(slice), .9 * max(slice), "FWHM: %f" % fwhm)

	# image display
	imagePlot.imshow(image, norm=lm(image.mean(), image.max(), clip='True'), cmap=cm.gray, origin="lower")
	imagePlot.set_xlim(xc - width, xc + width)
	imagePlot.set_ylim(yc - width, yc + width)

	# isophote plot
	if len(iso[0]) > 0:
		isophotePlot.scatter(*zip(*iso[0]), c="#ff00ff")
	if len(iso[1]) > 0:
		isophotePlot.scatter(*zip(*iso[1]), c="#0000ff")
	if len(iso[2]) > 0:
		isophotePlot.scatter(*zip(*iso[2]), c="#00ff00")
	if len(iso[3]) > 0:
		isophotePlot.scatter(*zip(*iso[3]), c="#ffff00")
	if len(iso[4]) > 0:
		isophotePlot.scatter(*zip(*iso[4]), c="#ff0000")

	isophotePlot.set_xlabel("X")
	isophotePlot.set_ylabel("Y")

	# aperture cdf
	apertureCDF.plot(radii, fluxes)

	apertureCDF.set_xlabel("Sample radius (px)")
	apertureCDF.set_ylabel("Flux (counts)")
	apertureCDF.set_title("Counts vs. Sample Radius")

	#pyplt.show()

	pyplt.tight_layout()

	pyplt.savefig("%s-analysis.pdf" % file)

if __name__ == '__main__':
	from sys import argv, exit

	try:
		fitsImage = fits.open(argv[1])
		image = fitsImage[0].data

		print("Finding center")
		(xc, yc, xp, yp, windowWidth, total) = center(image, float(argv[2]))

		radii = [x for x in range(0, int(windowWidth / 2))]

		fluxes = []

		print("Finding flux in radii")
		for radius in radii:
			fluxes.append(aperture(image, xc, yc, radius))

		print("Generating isophotes")
		iso = isophotes(image, float(argv[2]))

		print("Calculating FWHM")
		(slice, fwhm) = imageSlice(image, xc, yc, windowWidth)

		print("Generating plot")
		plot(image, xc, yc, windowWidth, slice, fwhm, iso, fluxes, radii, argv[1])

	except IOError:
		print("file does not exist")
		exit(-1)