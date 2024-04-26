[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11072638.svg)](https://doi.org/10.5281/zenodo.11072638)

# Comet 62P/Tsuchinshan
## A photometric investigation of an evolving dust production rate and other cometary properties at small distances post-perhelion

A collection of scripts and code written and used for my final year undergraduate research project within BSc. UCD Physics with Astronomy and Space Science.

Supervisor: Dr. Antonio Martin-Carillo

## Abstract

We presented photometric observations with Johnson-Cousins filters of comet 62P/Tsuchinshan 1
taken in March 2024. The dust rate proxy was determined for the standard reference aperture radius
$\rho = 10^4$ km, and was phase corrected to $0^\circ$ to yield an average R band $Af\rho$ = (73 ± 21) cm, B band
$Af\rho$ = (59 ± 5) cm, and V band $Af\rho$ = (70 ± 14) cm. We found absolute magnitudes in the B, V, and R
bands of 14.34 ± 0.16, 13.95 ± 0.23, and 13.59 ± 0.22 respectively, measured with a photometric aperture
radius $\rho = 10^4$ km. Dust colours were subsequently determined with colour indices B − V = 0.45 ± 0.3
and V − R = 0.3 ± 0.4 with reddening = (0.9 ± 2.7) % per 1000 Å between the B and V filters, and
reddening = (2 ± 4) % per 1000 Å between the V and R filters, from which the presence of gas emission
was inferred. An upper limit to the nucleus radius was determined through a coma correction method
and found to be $r_{\rm nucleus}$ ≤ 5.48 km.

_Keywords:_ Comets, 62P/Tsuchinshan 1, Photometry, Activity, Dust Colour, Nucleus Radius

## File Outline

Scripts are divided between two directories as follows:
- **DataAnalysisReview** contains Python scripts and notebooks created prior to the observations. These files were not used in the creation of figures for the final paper.
- **Scripts** contains Python scripts and notebooks to perform the data analysis and create the figures used in the final paper. These files in some places draw from ideas first formulated in **DataAnalysisReview**.

### ./Scripts/imageReduction.ipynb

Used to reduce raw telescope images using standard astronomy methods of bias/dark subtracting and flat field correcting.

### ./Scripts/imageStacking.ipynb

Used to stack reduced images from each day to create one coadded image centred on the comet.

### ./Scripts/activityAndColour.ipynb

Used to find the dust production rate of the comet, creating activity curves and activity evolution plots. Also used to find colour indices and reddening and related figures.

### ./Scripts/nucleusRadius.ipynb

Used to find an upper bound to the radius of the nucleus of the comet.

### ./Scripts/plotting.ipynb

Used to create other assorted plots such as: activity in the contex of other short and long period comets, and plots showing stacked images from different bands and nights together.

### ./Scripts/cometToolkit.py

A faux library containing common functions to avoid repetition across the files.

### Other

The other files **isophotes.ipynb** and **temporal.ipynb** were for investigating and creating plots which were ultimately omitted from the analysis.

## Data

All raw data and derived products can be accessed via Zenodo: 10.5281/zenodo.11072638
