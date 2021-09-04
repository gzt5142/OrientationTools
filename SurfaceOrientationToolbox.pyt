# -*- coding: utf-8 -*-
#
# License: GNU GPL v3 -- See 'LICENCE' file.
#
#

import numpy as np
import csv

from os.path import splitext, dirname, basename
from os import listdir
from datetime import date

from PIL import Image
import logging
from logging.handlers import RotatingFileHandler
import arcpy

from trf import surface, utils, shader, shadow
# GLOBALS:
# Enabling DEBUG will cause logging messages to be sent to a file handler as well as the arcpy message queue.
DEBUG = True


# Pre-amble -- used by all tools...
class ArcPyPassThru(logging.Handler):
    """
    Python Logging module handler to pass messages to the arcpy message queues. With this handler attached to
    a logger, you needn't call arcpy.AddMessage() separately -- nor should you if this is in use.
    """

    def __init__(self):
        super().__init__()

    def emit(self, record):
        msg = self.format(record)
        if record.levelno >= logging.ERROR:
            arcpy.AddError(msg)
            return
        if record.levelno >= logging.WARNING:
            arcpy.AddWarning(msg)
            return
        if record.levelno >= logging.INFO:
            arcpy.AddMessage(msg)


def initializeLoggers(a, lgr):
    """
    File-level logger instantiated.  this is done inside a function so as to keep the logger variable out of
    global scope.  If anybody wants to use this logger, I am compelling them to do this:
        logger = logging.getLogger(splitext(basename(__file__))[0])
    :return:
    """
    logger_ = logging.getLogger(splitext(basename(__file__))[0])
    if not logger_.handlers:
        logger_.addHandler(a)
        if DEBUG:
            logger_.addHandler(lgr)
            logger_.setLevel(logging.DEBUG)
        else:
            logger_.setLevel(logging.ERROR)
        logger_.debug(f"Toolbox Loaded: {date.today()}")


# Leaving these handlers in global scope so that various classes can use them
# to create their own loggers if they want.
logfileHandler = RotatingFileHandler(
    dirname(__file__) + "\\" + splitext(basename(__file__))[0] + ".log",
    maxBytes=1024 * 1024 * 2,  # 2MB log files
    backupCount=2
)
logfileHandler.setFormatter(
    logging.Formatter("%(levelname)-8s %(name)s %(asctime)s.%(msecs)03d %(message)s", datefmt='%H:%M:%S')
)

# Using the logger at INFO or above to pass the message at the
# appropriate level to the arcpy message queue.
# This will call arcpy.AddMessage() so you don't have to.
arcpyHandler = ArcPyPassThru()
arcpyHandler.setFormatter(
    logging.Formatter("%(asctime)s.%(msecs)03d >> %(message)s", datefmt='%H:%M:%S')
)

initializeLoggers(arcpyHandler, logfileHandler)



#  > > > > > > > > > > > > > > > Real Work Begins Here < < < < < < < < < < < < < < < <
class Toolbox(object):
    """
    Toolbox implementing terrain representation methods using Surface Normal Vectors as primary data model.
    """
    def __init__(self):
        self.label = "Surface Orientation Toolbox"
        self.description = self.__doc__
        self.alias = ''.join(filter(str.isalpha, splitext(basename(__file__))[0]))
        self.tools = [
            Study_DEM,
            Traditional_Hillshade,
            Traditional_Sky_Model,
            Soft_Hillshade,
            NLCD_Bump_Mapper,
            Prep_NLCD_Bumpmap_Mask
        ]


class ReliefTool(object):
    """
    This class does not implement anything. It only exists to sub-class .... this boilerplate tool class
    handles a lot of the code which would otherwise be replicated in every other tool below.  Note that
    this tool is not listed in the Toolbox.tools[] list, so it will not be seen by the user.
    """
    def __init__(self, **kwargs):
        fname = splitext(basename(__file__))[0]
        if 'toolname' in kwargs:
            self.label = kwargs.get('toolname').replace('_', ' ')
            self.LOG = logging.getLogger(fname + "." + kwargs.get('toolname'))
        else:
            self.label = __class__.__name__.replace('_', ' ')
            self.LOG = logging.getLogger(fname + "." + __class__.__name__)

        self.LOG.propagate = False
        if not self.LOG.handlers:
            self.LOG.addHandler(arcpyHandler)
            if DEBUG:
                self.LOG.addHandler(logfileHandler)
                self.LOG.setLevel(logging.DEBUG)
            else:
                self.LOG.setLevel(logging.INFO)
        # Anywhere in this tool/class, you can call self.LOG.<LEVEL>(message) and have the message
        # go to the right place: 1) arcpy msg queue and 2) if DEBUG, the logfile.

        self.description = self.__doc__
        #    ^^^^ We set this here to ensure that it has a value. It should be
        #         overridden in the subclass to be that class's doc string.


    def getParameterInfo(self):
        # These four parameters are fairly standard for any relief shading tool in this toolbox.
        # The subclassed tools can call super().getParaeterInfo() (i.e. this method) to get the boilerplate
        # parameter list then adjust.
        # Note the order, as it will be important later. The first line of the subclassed tool should be:
        # [inputDEM, lightAzimuth, lightElev, outputRaster] = super().getParameterInfo()

        inputDEM = arcpy.Parameter(
            displayName="Surface Raster (multiband)",
            name='input',
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input"
        )
        lightAzimuth = arcpy.Parameter(
            displayName="Light Azimuth",
            name="lightAz",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input"
        )
        lightAzimuth.value = 315

        lightElev = arcpy.Parameter(
            displayName="Light Elevation",
            name="lightEl",
            datatype="GPDouble",
            parameterType="Required",
            direction="Output"
        )
        lightElev.value = 45

        outputRaster = arcpy.Parameter(
            displayName="Output",
            name="output",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Output"
        )
        return [inputDEM, lightAzimuth, lightElev, outputRaster]

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        # We assume that the input multiband raster is the first parameter in the parameter list.
        # Maintain that convention.
        inputRaster = parameters[0].value
        parameters[0].clearMessage()
        if inputRaster:
            if not self.validSurface(inputRaster):
                parameters[0].setErrorMessage("This doesn't seem to be a 'studied' Surface...\nRun >>Study DEM<<")
            else:
                parameters[0].clearMessage()
        else:
            parameters[0].clearMessage()
        return

    def validSurface(self, path):
        """
        Runs from the tool validator functions.  Checks to see that the supplied path leads to one of our
        custom-generated multi-band rasters.  The multiband must have 4 bands, named DEM, Nx, Ny, and Nz.
        If not, returns false.
        :param path:
        :return:
        """
        oldEnv = arcpy.env.workspace
        arcpy.env.workspace = path
        bands = arcpy.ListRasters()
        arcpy.env.workspace = oldEnv

        if bands == ['DEM', 'Nx', 'Ny', 'Nz']:
            return True
        else:
            return False

    def readVectorArray(self, inputDEM):
        """
        Given a path to a multi-band, reads the Nx, Ny, and Nz named bands to return the 3D surface normal vectors.
        :param inputDEM: path to multiband
        :return: rank 3 numpy array
        """

        self.LOG.debug("Reading Vectors >> {}".format(basename(inputDEM)))
        arcpy.SetProgressor("step", "Reading from {}...".format(basename(inputDEM)), 0, 3, 1)
        try:
            Nx = arcpy.RasterToNumPyArray(inputDEM + r"\Nx")
            arcpy.SetProgressorPosition(1)
            Ny = arcpy.RasterToNumPyArray(inputDEM + r"\Ny")
            arcpy.SetProgressorPosition(2)
            Nz = arcpy.RasterToNumPyArray(inputDEM + r"\Nz")
            arcpy.SetProgressorPosition(3)
        except IOError as ex:
            self.LOG.exception("Could Not Read Surface Normal Data")
            raise ex
        S = np.stack([Nx, Ny, Nz], 0)
        return S

    def execute(self, parameters, messages):
        if DEBUG:
            arcpy.AddMessage("DEBUG logfile available in {}".format(basename(logfileHandler.baseFilename)))
        return

class Hello_World(ReliefTool):
    """
    Example tool to show how to subclass ReliefTool().
    """
    def __init__(self):
        super().__init__(toolname=__class__.__name__)
        self.description = self.__doc__
        self.category = "Primitives"

    def execute(self, parameters, messages):
        super().execute(parameters, messages)
        argv = {p.name: p for p in parameters}
        self.LOG.debug("Beginning Execution... ")   # should go to logfile only.
        self.LOG.info("Hello, World!")  # should go to logfile as well as the arcgis pro message window.
        self.LOG.debug("Done.") # logfile only.
        return

class Study_DEM(ReliefTool):
    """
    Processes a DEM to create a multi-band raster dataset representing the X,Y,Z components of the surface normal vectors.
    """
    def __init__(self):
        super().__init__(toolname=__class__.__name__)
        self.description = self.__doc__

    def getParameterInfo(self):
        # This is the one tool that doesn't use any of the usual parameters prototyped in ReliefTool.
        # Here, we just build our own from scratch.
        inputDEM = arcpy.Parameter(
            displayName="Input DEM",
            name='inDEM',
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input"
        )

        outputDataSet = arcpy.Parameter(
            displayName="Output",
            name="output",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output"
        )
        return [inputDEM, outputDataSet]


    def updateMessages(self, parameters):
        inputRaster = parameters[0].value
        parameters[0].clearMessage()
        if inputRaster:
            desc = arcpy.Describe(inputRaster)
            sr = desc.spatialReference
            if sr.type != "Projected":
                parameters[0].setErrorMessage("Input DEM must be in a Projected Coordinate System")
            else:
                parameters[0].clearMessage()
        else:
            parameters[0].clearMessage()
        return


    def execute(self, parameters, messages):
        super().execute(parameters, messages)
        argv = {p.name: p for p in parameters}

        self.LOG.debug("Reading DEM input: {}".format(argv['inDEM'].valueAsText))
        arcpy.SetProgressorLabel("Reading in DEM...")

        demInfo = arcpy.Describe(argv['inDEM'].valueAsText)
        tmpCorner = arcpy.Point(demInfo.extent.XMin, demInfo.extent.YMin)
        arcpy.env.outputCoordinateSystem = demInfo.spatialReference

        elevArray = arcpy.RasterToNumPyArray(arcpy.Raster( argv['inDEM'].valueAsText), nodata_to_value=0)

        self.LOG.debug("Calculating Surface Normals from Gradient.")
        arcpy.SetProgressorLabel("Calculating Gradient...")
        self.LOG.debug("Shape of elevation data: {}".format(elevArray.shape))

        sn = surface.normals_by_method(elevArray, demInfo.meanCellWidth, "N82")  # N82 method closely approximates the method Esri uses

        self.LOG.debug("Shape of surfacenormal output: {}".format(sn.shape))

        Nx = arcpy.NumPyArrayToRaster(sn[0,:,:], tmpCorner, demInfo.meanCellWidth, demInfo.meanCellHeight)
        Ny = arcpy.NumPyArrayToRaster(sn[1,:,:], tmpCorner, demInfo.meanCellWidth, demInfo.meanCellHeight)
        Nz = arcpy.NumPyArrayToRaster(sn[1,:,:], tmpCorner, demInfo.meanCellWidth, demInfo.meanCellHeight)
        elev = arcpy.NumPyArrayToRaster(elevArray, tmpCorner, demInfo.meanCellWidth, demInfo.meanCellHeight)

        self.LOG.debug("Composite to Multi-Band {}".format(basename(argv['output'].valueAsText)))
        arcpy.SetProgressorLabel("Composite to Multi-Band {}".format(argv['output'].valueAsText))
        arcpy.CompositeBands_management([elev, Nx, Ny, Nz], argv['output'].valueAsText)

        self.LOG.debug("Renaming to DEM, Nx, Ny, Nz")
        arcpy.SetProgressor("step", "Renaming Bands...", 0, 4, 1)
        arcpy.Rename_management(argv['output'].valueAsText + "\\Band_1", "DEM")
        arcpy.SetProgressorPosition(1)
        arcpy.Rename_management(argv['output'].valueAsText + "\\Band_2", "Nx")
        arcpy.SetProgressorPosition(2)
        arcpy.Rename_management(argv['output'].valueAsText + "\\Band_3", "Ny")
        arcpy.SetProgressorPosition(3)
        arcpy.Rename_management(argv['output'].valueAsText + "\\Band_4", "Nz")
        arcpy.SetProgressorPosition(4)

        self.LOG.debug("Finished.")
        return

class Traditional_Hillshade(ReliefTool):
    """
    Clamps output of the Lambertian shader to [0..255], just as the traditional hillshade from Esri, Grass, GDAL, etc.
    Unlike pure Lambert, this output is unsigned bytes (uint8) to match what others produce.
    """
    def __init__(self):
        super().__init__(toolname=__class__.__name__)
        self.description = self.__doc__
        self.category = "Shaded Relief"

    def getParameterInfo(self):
        [inputDEM, lightAzimuth, lightElev, outputRaster] = super().getParameterInfo()
        modelShadows = arcpy.Parameter(
            displayName = "Model Shadows?",
            name='shadows',
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input"
        )
        modelShadows.value = False
        return [inputDEM, lightAzimuth, lightElev, modelShadows, outputRaster]

    def execute(self, parameters, messages):
        super().execute(parameters, messages)
        argv = {p.name: p for p in parameters} # convenience dict

        ## Preamble...
        inputDEM = argv['input'].valueAsText
        arcpy.SetProgressorLabel("Setting Up...")
        self.LOG.debug("Read Inputs...")
        demInfo = arcpy.Describe(inputDEM + r"\DEM")
        tmpCorner = arcpy.Point(demInfo.extent.XMin, demInfo.extent.YMin)
        arcpy.env.outputCoordinateSystem = demInfo.spatialReference
        surfaceNormals = self.readVectorArray(inputDEM)

        # Ready to do work
        arcpy.SetProgressor("default", "Calculating Shade Values...")
        self.LOG.debug("Shading...")

        L = utils.lightVector(argv['lightAz'].value, argv['lightEl'].value)
        hs = shader.lambert(surfaceNormals, L)
        self.LOG.debug(f"Convert {hs.dtype} to uint8 ...")
        ## All negative values set to zero...
        hs[hs < 0] = 0
        # Scale to range [0-255]  -- but we are still a floating point data type!
        hs *= 255.0
        # Traditional hillshades are typically returned as unsigned 8-bit chars/bytes. Must cast
        # to uint8 at some point before saving....
        clamped = hs.astype('uint8')
        self.LOG.debug(f"Output Clamped to [0..255] Min={clamped.min()} / Max={clamped.max()}")

        if argv['shadows'].value == True:
            arcpy.SetProgressorLabel("Modeling Shadows...")
            self.LOG.debug(f"Shadows ...")
            elevArray = arcpy.RasterToNumPyArray(inputDEM + r"\DEM")
            shadowArray = shadow.shadowLine(elevArray,
                                                 argv['lightAz'].value,
                                                 argv['lightEl'].value,
                                                 demInfo.meanCellWidth
                                                 )
            clamped[shadowArray > 0] = 0
        arcpy.SetProgressor("default", "Writing Output to {} ...".format(basename(argv['output'].valueAsText)))
        self.LOG.debug("Writing Outputs...")
        outRaster = arcpy.NumPyArrayToRaster(clamped, tmpCorner, demInfo.meanCellWidth, demInfo.meanCellHeight)
        outRaster.save(argv['output'].valueAsText)
        return

class Soft_Hillshade(ReliefTool):
    """
    Produces shaded output with BVs in range 0 to 1.  values below 0.5 are self-shaded.
    This is essentially the output from lambert, re-scaled to fit into the range [0-1] rather than the
    cosine output of [-1 - 0].
    """
    def __init__(self):
        super().__init__(toolname=__class__.__name__)
        self.description = self.__doc__
        self.category = "Shaded Relief"

    def getParameterInfo(self):
        [inputDEM, lightAzimuth, lightElev, outputRaster] = super().getParameterInfo()
        modelShadows = arcpy.Parameter(
            displayName = "Model Shadows?",
            name='shadows',
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input"
        )
        modelShadows.value = False
        return [inputDEM, lightAzimuth, lightElev, modelShadows, outputRaster]

    def execute(self, parameters, messages):
        super().execute(parameters, messages)
        argv = {p.name: p for p in parameters} # convenience dict

        # Preamble...
        inputDEM = argv['input'].valueAsText
        arcpy.SetProgressorLabel("Setting Up...")
        self.LOG.debug("Read Inputs...")
        demInfo = arcpy.Describe(inputDEM + r"\DEM")
        tmpCorner = arcpy.Point(demInfo.extent.XMin, demInfo.extent.YMin)
        arcpy.env.outputCoordinateSystem = demInfo.spatialReference
        surfaceNormals = self.readVectorArray(inputDEM)

        # Ready to do work
        arcpy.SetProgressor("default", "Calculating Shade Values...")
        self.LOG.debug("Shading...")

        L = utils.lightVector(argv['lightAz'].value, argv['lightEl'].value)
        hs = shader.lambert(surfaceNormals, L)  # lambert output is -1 to 1...
        hs =  (hs + 1 ) / 2  # rescale to fit 0 to 1

        if argv['shadows'].value == True:
            self.LOG.debug("Casting Shadows...")
            arcpy.SetProgressorLabel("Casting Shadows...")
            elevArray = arcpy.RasterToNumPyArray(inputDEM + r"\DEM")
            shadowArray = shadow.shadowLine(elevArray,
                                                 argv['lightAz'].value,
                                                 argv['lightEl'].value,
                                                 demInfo.meanCellWidth
                                                 )

            hs[shadowArray > 0] *= 0.5

        arcpy.SetProgressor("default", "Writing Output to {} ...".format(basename(argv['output'].valueAsText)))
        self.LOG.debug("Writing Outputs...")
        outRaster = arcpy.NumPyArrayToRaster(hs, tmpCorner, demInfo.meanCellWidth, demInfo.meanCellHeight)
        outRaster.save(argv['output'].valueAsText)
        return

class Sky_Tool(ReliefTool):
    """
    Empty class.  Sky model tools to subclass this.   Used mostly to put routines in one place
    rather than duplicating code.
    """
    def __init__(self, **kwargs):
        if 'toolname' in kwargs:
            super().__init__(toolname=kwargs.get('toolname'))
        else:
            super().__init__(toolname=__class__.__name__)
        self.description = self.__doc__
        self.category = "Shaded Relief"

    def getParameterInfo(self):
        [inputDEM, _, _, outputRaster] = super().getParameterInfo()
        skyviewConfigFile = arcpy.Parameter(
            displayName="Sky Configuration File:",
            name='config',
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        )
        try:
            dname = dirname(__file__) + "\\SkyConfigFiles\\"
            self.LOG.debug("Looking for sky configs in {}... ".format(dname))
            skyviewConfigFile.filter.list = sorted([f for f in listdir(dname) if f.endswith('.txt')])
            self.LOG.debug("...found {}: {}".format(len(skyviewConfigFile.filter.list), skyviewConfigFile.filter.list))
        except:
            raise

        shadowLength = arcpy.Parameter(
            displayName="Shadow-Length",
            name='ShadowLength',
            datatype="GPDouble",
            parameterType="Required",
            direction="Input"
        )
        shadowLength.value = 5
        shadowLength.filter.list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        shadowOnly = arcpy.Parameter(
            displayName="Shade/Shadow Mode:",
            name='ShadeShadow',
            datatype="GPString",
            parameterType="Optional",
            direction="Input"
        )
        shadowOnly.value = "Shade & Shadow"
        shadowOnly.filter.list = ["Shade & Shadow", "Shade Only", "Shadow Only"]

        return [inputDEM, skyviewConfigFile, shadowLength, shadowOnly, outputRaster]

    def sky_light_list(self, f):
        self.LOG.debug(f"Reading Sky Config File: {f}")
        lts = []
        try:
            with open(f, "r") as skyConfigFile:
                csvReadFile = csv.reader(skyConfigFile)
                for line in csvReadFile:
                    # brute force -- if this line does not have three comma-separated values, we skip it.
                    if (len(line) == 0) or (len(line) < 3):
                        continue
                    if line[0].startswith("Format"):  # special case... a comment line in the header has 2 commas in it.
                        continue
                    lts.append(line)
        except IOError as e:
            self.LOG.exception("I/O Error: {}".format(e))
            raise e
        return lts

    def aggregateShadeLayers(self, DEMarray, demInfo, N, lightList, shade_shadow, clampOutput):
        outputArray = np.zeros(DEMarray.shape)
        arcpy.SetProgressor("step", f"Processing {len(lightList)} Illumination Sources...", 0, len(lightList), 1)
        for (i, (azimuth, elev, wt)) in enumerate(lightList):
            self.LOG.debug(f"...Light {i}: az={azimuth}, el={elev}, wt={wt}")
            L = utils.lightVector(float(azimuth), float(elev))
            arcpy.SetProgressorLabel(f"Light {i} of {len(lightList)}:  Azimuth={azimuth}, Elevation={elev}")
            arcpy.SetProgressorPosition(i)

            castshadow = True
            shade = True
            if (shade_shadow == "Shade & Shadow"):
                shade = True
                castshadow = True
            if (shade_shadow == "Shadow Only"):
                shade = False
                castshadow = True
            if (shade_shadow == "Shade Only"):
                shade = True
                castshadow = False

            if (castshadow):
                shadowArray = shadow.shadowLine(DEMarray, float(azimuth), float(elev),  demInfo.meanCellWidth )
            else:
                shadowArray = np.zeros(DEMarray.shape)

            if (shade):
                hs = shader.lambert(N, L)
            else:
                hs = np.ones(DEMarray.shape)


            if (clampOutput):
                hs[hs<=0] = 0
                hs[shadowArray > 0] = 0
                hs *= 255.0
                outputArray = outputArray + (hs * float(wt))

            else:
                hs = (hs + 1) / 2
                hs[shadowArray > 0] *= 0.5
                outputArray = outputArray + (hs * float(wt))
        self.LOG.debug(f"Processed {i+1} light vectors.")
        return outputArray

class Traditional_Sky_Model(Sky_Tool):
    """
    Kennelly, P. J., & Stewart, A. J. (2014). General sky models for illuminating terrains. International Journal
    of Geographical Information Science, 28(2), 383–406. https://doi.org/10.1080/13658816.2013.848985
    """
    def __init__(self, **kwargs):
        super().__init__(toolname=__class__.__name__)
        self.description = self.__doc__
        self.category = "Shaded Relief"


    def execute(self, parameters, messages):
        argv = {p.name: p for p in parameters}
        inputDEM = argv['input'].valueAsText
        arcpy.SetProgressorLabel("Setting Up...")
        self.LOG.debug("Read Inputs...")
        demInfo = arcpy.Describe(inputDEM + r"\DEM")
        tmpCorner = arcpy.Point(demInfo.extent.XMin, demInfo.extent.YMin)
        arcpy.env.outputCoordinateSystem = demInfo.spatialReference
        elevArray = arcpy.RasterToNumPyArray(inputDEM + r"\DEM").astype('float32') * argv['ShadowLength'].value

        surfaceNormals = self.readVectorArray(inputDEM)

        confFile = dirname(__file__) + "\\SkyConfigFiles\\" + argv['config'].valueAsText
        lights = self.sky_light_list(confFile)

        self.LOG.debug("Shade_Shadow = '{}'".format(argv['ShadeShadow'].value))
        outputArray = self.aggregateShadeLayers(elevArray,
                                                demInfo,
                                                surfaceNormals,
                                                lights,
                                                argv['ShadeShadow'].value,
                                                True)

        outRaster = arcpy.NumPyArrayToRaster(outputArray.astype('uint8'),
                                             tmpCorner,
                                             demInfo.meanCellWidth,
                                             demInfo.meanCellHeight
                                             )
        arcpy.SetProgressor("default", "Writing Output to {} ...".format(basename(argv['output'].valueAsText)))
        outRaster.save(argv['output'].valueAsText)
        return

class Bump_Tool(ReliefTool):
    def __init__(self, **kwargs):
        if 'toolname' in kwargs:
            super().__init__(toolname=kwargs.get('toolname'))
        else:
            super().__init__(toolname=__class__.__name__)
        self.description = self.__doc__

    def reScale(self, bumpMap):
        logger = logging.getLogger(splitext(basename(__file__))[0])
        logger.debug("...Casting {} to range -1 to 1".format(bumpMap.dtype))
        if bumpMap.dtype == 'uint8':
            return ((bumpMap / (2 ** 8 - 1)) * 2) - 1
        if bumpMap.dtype == 'uint16':
            return ((bumpMap / (2 ** 16 - 1)) * 2) - 1
        if bumpMap.dtype == 'float32':
            return (bumpMap - 0.5) * 2
        if bumpMap.dtype == 'float64':
            return (bumpMap - 0.5) * 2
        return bumpMap

    def consolidateBumpMaps(self, extentShape, nlcdArray, t):
        """
        :param extentShape: The tuple representing the shape of the output array.  Should be (band,row,col), with three bands.
        :param nlcdArray:
        :param t:
        :return:
        """
        logger = logging.getLogger(splitext(basename(__file__))[0])
        logger.debug("Consolidating bump maps...")

        bmMaster = np.zeros(extentShape)
        bmMaster[2, :,:] = 1.0 # all vectors in the bump-map master are now [0,0,1]
        logger.debug("bmMaster shape is {}".format(bmMaster.shape))
        logger.debug("NLCD shape is {}".format(nlcdArray.shape))

        pathPrefix = dirname(__file__) + "\\BumpMaps\\"
        arcpy.SetProgressor("step", "Tiling Bump Map:", 0, len(t)+1, 1)
        for (count, row) in enumerate(t):
            arcpy.SetProgressorPosition(count)
            maskValue = row[0]
            bumpmapFile = pathPrefix + row[1]
            logger.debug("Processing bump map: {} of {}: {}".format(count, len(t), bumpmapFile))
            image = Image.open(bumpmapFile)
            Bm = np.asarray(image)
            # IMPORTANT NOTE:  An RGB image is read in to an array of shape (row, column, band).
            # In manipulating multi-band rasters in arcpy, the band is the first index.  Needs to be
            # re-shaped to (band, row, col).  We'll do that after tiling the image...

            logger.debug(" ...Tiling {}x{} Image: {}".format(Bm.shape[0], Bm.shape[1], row[1]))
            _, x, y = bmMaster.shape
            x_tiles = (x // Bm.shape[0]) + 2
            y_tiles = (y // Bm.shape[1]) + 2
            img = np.tile(Bm, (x_tiles, y_tiles, 1))
            center_x = img.shape[0] // 2
            center_y = img.shape[1] // 2
            Bm = img[center_x - (x // 2):center_x + (x // 2) + 1, center_y - (y // 2):center_y + (y // 2) + 1, :]
            Bm = Bm[0:x, 0:y, 0:3] # discard alpha channel, if present
            B = self.reScale(Bm)
            B = np.moveaxis(B, 2, 0) ##<<<<< change from (row,col,band) to (band,row,col), to match bmMaster.
            logger.debug("B shape is {}".format(B.shape))
            (x,y) = nlcdArray.shape
            n = np.broadcast_to(nlcdArray, (3, x, y))
            bmMaster[n == maskValue] = B[n == maskValue]

        logger.debug("Master Bump Map Assembled... ")
        return  bmMaster

    def applyBumpMap(self, S, B):
        """
        :param S: Surface normal vector array.  (band, row, column)
        :param B: Bump map surface normal vectors (band, row, column)
        :return: Bumpified surface normal vectors for surface
        """
        logger = logging.getLogger(splitext(basename(__file__))[0])
        logger.debug("Bump Map >> S shape is {}, B shape is {}".format(S.shape, B.shape))
        Uz = S[0,:, :] / S[2, :, :]
        mag = np.sqrt(1 + Uz ** 2)
        U = np.stack([1 / mag, np.zeros((S.shape[1], S.shape[2])), -Uz / mag], 0)

        Vz = S[1, :, :] / S[2, :, :]
        mag = np.sqrt(1 + Vz ** 2)
        V = np.stack([np.zeros((S.shape[1], S.shape[2])), 1 / mag, -Vz / mag], 0)

        # N is a straightforward change-of-basis calculation.  We are going from tangent space
        # to world space.  Tangent space is defined by U, V, and S vectors:
        #    - U is i-hat (vector representing one unit in "x" direction;
        #    - V is j-hat (vector representing one unit in "y" direction;
        #    - S is the original surface normal, a.k.a. k-hat -- one unit in the "z" direction.
        # Dot each of these vectors with the bump-map vector (which is in tangent space)
        # to get the i,j,k components of that same vector in world space:
        Nx = (U[0, :, :] * B[0, :, :]) + (V[0,:, :] * B[1,:, :]) + (S[0, :, :] * B[2, :, :])
        Ny = (U[1, :, :] * B[0, :, :]) + (V[1,:, :] * B[1,:, :]) + (S[1, :, :] * B[2, :, :])
        Nz = (U[2, :, :] * B[0, :, :]) + (V[2,:, :] * B[1,:, :]) + (S[2, :, :] * B[2, :, :])

        mag = np.sqrt(Nx ** 2 + Ny ** 2 + Nz ** 2)  # In theory, N is already normalized.  This guarantees it.
        # N is the surface normal vector S, with bump map B applied to it.
        logger.debug("Bumped.")
        return np.stack([Nx / mag, Ny / mag, Nz / mag], 0)


class NLCD_Bump_Mapper(Bump_Tool):
    def __init__(self, **kwargs):
        super().__init__(toolname=__class__.__name__)
        self.description = self.__doc__
        self.category = "BumpMaps"

    def getParameterInfo(self):
        [inputDEM, lightAzimuth, lightElev, outputRaster] = super().getParameterInfo()
        nlcd = arcpy.Parameter(
            displayName="NLCD:",
            name='NLCD',
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input"
        )
        vt = arcpy.Parameter(
            name='lc_table',
            displayName='Land Cover Table',
            datatype='GPValueTable',
            direction='Input',
            parameterType='Optional'
        )
        vt.parameterDependencies = [inputDEM.name, nlcd.name]
        vt.columns = [
            ['GPLong', 'NLCD Code'],
            ['GPString', 'Bumpmap Tile'],
        ]
        vt.enabled = False
        return [inputDEM, lightAzimuth, lightElev, outputRaster, nlcd, vt]

    def updateMessages(self, parameters):
        super().updateMessages(parameters)
        if parameters[4].valueAsText:
            nlcd = arcpy.Raster(parameters[4].valueAsText)
            parameters[4].clearMessage()
            if not nlcd.isInteger:
                parameters[4].setErrorMessage("Land Use layer must be integer pixel type.")
                self.LOG.debug(f"REJECTED: Tried to set land use layer with a '{nlcd.pixelType}' pixel type")
            else:
                parameters[4].clearMessage()
        return

    def updateParameters(self, params):
        if params[4].valueAsText:
            if params[5].enabled:
                pass
            else:
                defaultValues = [11, 12, 22, 23, 24, 31, 41, 42, 43, 51, 52, 71, 72, 73, 74, 81, 82, 90, 95]
                self.LOG.debug("Enabling NLCD value table... ")
                params[5].enabled = True
                self.LOG.debug("Scanning {} for unique land use codes...".format(params[4].valueAsText))
                uniqValues = getUniqueValues(params[4].valueAsText)
                if len(uniqValues) > 0:
                    self.LOG.debug(f"Found {len(uniqValues)} NLCD codes.")
                    params[5].filters[0].list = uniqValues
                else:
                    self.LOG.debug(f"Found none.  Setting to default list.")
                    params[5].filters[0].list = defaultValues
                params[5].filters[1].type = "ValueList"

                try:
                    dname = dirname(__file__) + "\\BumpMaps\\"
                    self.LOG.debug("Looking for bump maps in {}... ".format(dname))
                    params[5].filters[1].list = sorted([f for f in listdir(dname) if f.endswith('.tiff')] )
                    self.LOG.debug("...found {}".format(len(params[5].filters[1].list)))
                except:
                    raise
        else:
            params[5].enabled = False

        return

    def execute(self, parameters, messages):
        argv = {p.name: p for p in parameters}
        inputDEM = argv['input'].valueAsText
        arcpy.SetProgressorLabel("Setting Up...")
        self.LOG.debug("Read Inputs...")
        demInfo = arcpy.Describe(inputDEM + r"\DEM")
        tmpCorner = arcpy.Point(demInfo.extent.XMin, demInfo.extent.YMin)
        arcpy.env.outputCoordinateSystem = demInfo.spatialReference

        surfaceNormals = self.readVectorArray(inputDEM)

        nlcd = arcpy.RasterToNumPyArray(argv['NLCD'].valueAsText)
        nlcdInfo = arcpy.Describe(argv['NLCD'].valueAsText)

        ## nlcd raster must have same geometry and projection as the surface dataset. . .
        if demInfo.spatialReference.name != nlcdInfo.spatialReference.name:
            self.LOG.error("Spatial Reference Does Not Match")

        L = utils.lightVector(argv['lightAz'].value, argv['lightEl'].value)

        xDim = min(surfaceNormals.shape[1], nlcd.shape[0])
        yDim = min(surfaceNormals.shape[2], nlcd.shape[1])
        nlcd = nlcd[0:xDim, 0:yDim]
        self.LOG.debug("NLCD read from {}, with shape {}".format(argv['NLCD'].valueAsText, nlcd.shape))

        bmMaster = self.consolidateBumpMaps(surfaceNormals.shape, nlcd, argv['lc_table'].values)

        arcpy.SetProgressor("default", "Applying Assembled Bump Map ...")
        N = self.applyBumpMap(surfaceNormals, bmMaster)

        self.LOG.debug(" ...Shading Bumpmapped normal field")
        bmHS = shader.lambert(N, L)  # output is -1 to 1
        bmHS = (bmHS + 1) / 2 # scale to fit 0-1 range
        bmHS = bmHS[0:xDim,0:yDim]
        self.LOG.debug("DONE")

        arcpy.SetProgressor("default", "Writing Output to {} ...".format(basename(argv['output'].valueAsText)))
        hs = arcpy.NumPyArrayToRaster(bmHS, tmpCorner, demInfo.meanCellWidth, demInfo.meanCellHeight)
        hs.save(argv['output'].valueAsText)
        return

class Prep_NLCD_Bumpmap_Mask(ReliefTool):
    def __init__(self):
        super().__init__(toolname=__class__.__name__)
        self.description = self.__doc__
        self.category = "BumpMaps"

    def getParameterInfo(self):
        [inputDEM, _, _, outputRaster] = super().getParameterInfo()
        nlcd = arcpy.Parameter(
            displayName="NLCD:",
            name='NLCD',
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input"
        )
        return [inputDEM, nlcd, outputRaster]

    def execute(self, parameters, messages):
        argv = {p.name: p for p in parameters}
        inputDEM = argv['input'].valueAsText

        arcpy.SetProgressorLabel("Setting Up...")
        self.LOG.debug("Read Inputs...")

        demInfo = arcpy.Describe(inputDEM + r"\DEM")
        arcpy.env.outputCoordinateSystem = demInfo.spatialReference
        arcpy.env.extent = demInfo.extent
        arcpy.env.snapRaster = inputDEM + "r\DEM"
        envelope = str(demInfo.extent.XMin) + " " + str(demInfo.extent.YMin) + " " + str(demInfo.extent.XMax) + " " + str(demInfo.extent.YMax)

        self.LOG.debug("Resampling...")
        arcpy.Resample_management(argv['NLCD'].valueAsText, "memory\\tmp_NLCD", demInfo.meanCellWidth, "NEAREST")
        self.LOG.debug("Clipping...")
        arcpy.Clip_management("memory\\tmp_NLCD", envelope, argv['output'].valueAsText, "#", "#", "#", "MAINTAIN_EXTENT")
        self.LOG.debug("Cleaning Up...")
        arcpy.Delete_management("memory\\tmp_NLCD")
        return

#  > > > > > > > > > > > > > > > Useful utility functions < < < < < < < < < < < < < < < <

def getUniqueValues(P):
    r = arcpy.Raster(P)
    if r.isInteger:
        a = arcpy.RasterToNumPyArray(r)
        return np.unique(a.ravel()).tolist()
    else:
        return []
