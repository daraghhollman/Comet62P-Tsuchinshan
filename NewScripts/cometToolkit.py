# Importing packages

from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import mad_std

from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, aperture_photometry
from photutils.background import Background2D

from natsort import natsorted

import matplotlib.pyplot as plt

import numpy as np


def GetImage(path):
    """
    Loads the image from a fits file
    """

    return fits.open(path)[0].data


def SaveFits(data, outputPath, header=None):
    """
    Saving fits files
    """

    newFits = fits.PrimaryHDU(data=data)
    #newFits.data = data.astype(int)

    print(np.max(data))
    print(np.max(newFits.data))

    if header is not None:
        newFits.header = header

    newFits.writeto(outputPath, overwrite=True)
    print(f"FITS image saved to: {outputPath}")



def PlotFits(path, wcsPath="", vmin=None, vmax=None,
             cmap="binary"):
    
    """
    Plots a fits image file

    Positional Arguments:
    path -- the .fits path to load from

    Keyword arguments:
    wcsPath -- path to wcs fits file
    vmin -- colourbar minimum pixel counts
    vmax -- colourbar maxiumum pixel counts
    cmap -- matplotlib cmap
    """


    # Extracting image from fits
    image = GetImage(path)

    # If wcs is availble, use wcs projection
    
    if wcsPath != "":
        wcs = WCS(fits.open(wcsPath)[0].header)

        plt.subplot(projection=wcs)

    
    # Plotting image with colourbar
    plt.imshow(image, vmin=vmin, vmax=vmax, cmap=cmap)
    plt.colorbar(label="Counts")

    # Setting axis labels
    if wcsPath != "":
        plt.xlabel("Right Ascension")
        plt.ylabel("Declination")

    else:
        plt.xlabel("CCD X")
        plt.ylabel("CCD Y")


##### IMAGE MANIPULATION #####

def RollImage(image, rollTuple):
    return np.roll(image, rollTuple, axis=(1,0))



###### IMAGE REDUCTION #######

def ReduceImage(inputPath, outputPath, biasPath, darkPath, flatPath, showPlot=False):

    # Object image minus bias minus darks (bias subtracted) and divided by the master flats

    fitsFile = fits.open(inputPath)
    biasImage = fits.open(biasPath)[0].data
    darkImage = fits.open(darkPath)[0].data
    flatImage = fits.open(flatPath)[0].data

    objectImage = fitsFile[0].data
    header = fitsFile[0].header

    reducedImage = (objectImage - biasImage - darkImage) / flatImage

    #SaveFits(reducedImage, outputPath, header=header)

    newFits = fits.PrimaryHDU(data=reducedImage)
    newFits.writeto(outputPath, overwrite=True)

    if showPlot:
        plt.imshow(reducedImage, vmax=16000)

    return


###### STACKING ########

def SearchStars(image, fwhm=8., threshold=4, showPlot=False):

    medianAbsoluteDeviation = mad_std(image)
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold*medianAbsoluteDeviation)
    sources = daofind(image)

    if (showPlot):
        plt.imshow(image)
        plt.scatter(sources['xcentroid'], sources['ycentroid'], alpha=0.5, color="orange")
        plt.xlabel("Arbitrary X")
        plt.ylabel("Arbitrary Y")
        plt.gca().invert_yaxis()

    return sources




def FindCometCentre(path, filter, maxCentreDistance=500,
                    showPlot=False):
    
    image = GetImage(path)

    # FInd all sources in sky
    sources = SearchStars(image, fwhm=8, threshold=16, showPlot=False)

    # Remove sources based on peak counts and distance from centre

    # Values calibrated manually per filter
    match filter:
        case "V":
            maxPeak = 1400
            minPeak = 700

        case "R":
            maxPeak = 2400
            minPeak = 2000

        case "B":
            maxPeak = 1000
            minPeak = 200


    removeIndices = []

    for i in range(len(sources)):

        # remove based on peak
        if sources[i]["peak"] > maxPeak or sources[i]["peak"] < minPeak:
            removeIndices.append(i)

        # remove based on position
        if np.sqrt( (sources[i]["xcentroid"] - 1024)**2 + (sources[i]["ycentroid"] - 1024)**2 ) > maxCentreDistance:
            removeIndices.append(i)

    updatedSources = sources
    updatedSources.remove_rows(removeIndices)

    
    if showPlot:
        
        match filter:
            case "V":
                upperThreshold = 1200
                lowerThreshold = 200

            case "R":
                upperThreshold = 1000
                lowerThreshold = 500

            case "B":
                upperThreshold = 1000
                lowerThreshold = 200

        thresholdIndices = np.where((image > lowerThreshold)  & (image < upperThreshold), 1, 0)

        plt.imshow(thresholdIndices)

        for source in sources:
            plt.scatter(source["xcentroid"], source["ycentroid"])


    # Add up pixels around the source
    # The comet will have brighter surrounding pixels for sources
    # of equal counts
    r = 50
    brightestCounts = 0
    cometId = None

    # Add pixels around source in radius
    for i, source in enumerate(updatedSources):

        # Loop through and add pixels in box size r
        x = round(source["xcentroid"])
        y = round(source["ycentroid"])

        counts = np.sum(image[x - r : x + r, y - r : y + r])

        if counts > brightestCounts:
            brightestCounts = counts
            cometId = i

    
    cometCentre = (updatedSources[i]["xcentroid"], updatedSources[i]["ycentroid"])

    return cometCentre




def MedianStack(imagesList, preserveAvgCounts=False):
    """
    Correctly median stacks a list of images by multiplying by the sqrt of the number of images
    """

    if preserveAvgCounts:
        return np.median(imagesList, axis=0)

    else:
        return np.median(imagesList, axis=0) * np.sqrt(len(imagesList))


def StackImages(paths):

    imagesToStack = []

    for path in paths:

        image = fits.open(path)[0].data

        imagesToStack.append(image)

    stackedImage = MedianStack(imagesToStack, preserveAvgCounts=True)

    return stackedImage


##### APERTURE PHOTOMETRY #####

def DetermineStarZeroPoint(stackedImage, numStacks, coordinatesList, calibratedMagnitudes, apertureRadius=30, tolerance=20, showPlot=False):
    
    """
    For a list of star coordinates, and corresponding list
    of calibrated magnitudes, returns the zeropoint for each star.

    POSITIONAL ARGUMENTS:
    stackedImage -- fits image data of stacked images with respect to background (must be in comet images)
    numStacks -- number of images used to create the stacked image
    coordinatesList -- a list of tuples with star x and y in pixel coordinates
    calibratedMagnitudes -- a list of calibration magnitudes corresponding to the stars

    KEYWORD ARGUMENTS:
    apertureRadius -- The DAOStarFinder aperture radius to do calculations in
    tolerance -- The pixel tolerance on the coordinates list
    showPlot -- Plots the sources from DAOStarFinder with all stars within tolerance
    """

    sources = SearchStars(stackedImage, showPlot=showPlot)

    zeroPoints = []

    # Finding stars
    for i, starCoords in enumerate(coordinatesList):
        starX = starCoords[0]
        starY = starCoords[1]

        referenceStarIndices = np.where((abs(sources["xcentroid"] - starX) < tolerance) & (abs(sources["ycentroid"] - starY) < tolerance))

        if showPlot:        
            plt.scatter(sources[referenceStarIndices]["xcentroid"], sources[referenceStarIndices]["ycentroid"], alpha=0.5, color="orange")

        #print(sources[referenceStarIndices])
        if len(sources[referenceStarIndices]) > 1:
            raise Exception("Too many stars found within tolerance for star: " + str(i))

        if len(sources[referenceStarIndices]) == 0:
            raise Exception("No stars found within tolerance for star: " + str(i))
        
        # Create aperture
        aperture = CircularAperture((sources[referenceStarIndices][0]["xcentroid"], sources[referenceStarIndices][0]["ycentroid"]), r=apertureRadius)

        # Find background
        background = Background2D(stackedImage, 50).background

        phot_table = aperture_photometry(stackedImage - background, aperture)
        
        observedMagnitude = -2.5 * np.log10( phot_table["aperture_sum"] / (120 * numStacks))

        
        zeroPoint = calibratedMagnitudes[i] - observedMagnitude
        zeroPoints.append(zeroPoint[0])

    return zeroPoints




##### ACTIVITY #####

def PlotAfrho(apertureRange, activityValuesV, activityValuesR, activityValuesB, colours=["green", "red", "blue"],
                labels=["V", "R", "B"]):
    """
    Plots activity curve as a function of aperture
    
    Positional Arguments:
    apertureRange -- list of apertures
    activityValues -- list of corresponding activities
    """

    fig, ax = plt.subplots()

    for i, activityValues in enumerate([activityValuesV, activityValuesR, activityValuesB]):
        ax.plot(apertureRange, activityValues, color=colours[i], label=labels[i])

    ax.set_xlabel(r"$\rho$ [km]")
    ax.set_ylabel(r"$Af\rho$ [cm]")
    ax.legend()

