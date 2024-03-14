# Importing packages

from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import mad_std

from photutils.detection import DAOStarFinder

from natsort import natsorted

import matplotlib.pyplot as plt

import numpy as np


def GetImage(path):
    """
    Loads the image from a fits file
    """

    return fits.open(path)[0].data


def SaveFits(data, outputPath):
    """
    Saving fits files
    """

    newFits = fits.PrimaryHDU()
    newFits.data = data

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

    
    # Flip declination axis
    plt.gca().invert_yaxis()




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
            maxPeak = 2400
            minPeak = 1800

        case "R":
            maxPeak = 3200
            minPeak = 2000

        case "B":
            maxPeak = 1350
            minPeak = 1200


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
                upperThreshold = 1600
                lowerThreshold = 1150

            case "R":
                upperThreshold = 1200
                lowerThreshold = 1040

            case "B":
                upperThreshold = 1200
                lowerThreshold = 1040

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
