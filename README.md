## Orientation Tools


### Python Toolbox for ArcGIS Pro: Surface Orientation Toolbox

NOTE: This toolbox requires the Pillow image processing package.  Pillow is included in the ArcGIS Pro distribution as of
version 2.7.  If you are using an earlier version of Pro, you will need to clone your python environment and install the
Pillow package manually. 


### Toolbox Contents

* __Study DEM__ -- Tool to process DEM inputs and generate the multi-band raster dataset representing the surface normal 
  vectors used by all of the tools listed below.  Surface normal vectors are derived from gradient obtained via 
  a Sobel operator.  The output from this tool should be placed in a geodatabase.  It will be a 4-band raster,
  with the bands being:
    * DEM -- A copy of the elevation data
    * Nx -- The X component of the surface normal vector
    * Ny -- The Y component of the surface normal vector
    * Nz -- The Z component of the surface normal vector

  The output from `Study DEM` will be input for the following tools: 
  
* _<u>Shaded Relief</u>_ ToolSet
    * __Hard Hillshade__ -- Produces a hillshade layer similar to that produced by other major players (Esri,
      GDAL, GRASS, etc). Output is clampted to unsigned integers from the range 0 to 255. 
      All negative values from the Lambertian
      shader are reset to zero for the traditional output. Option for shadows to be modeled. Any cell found to be in
      shadow is set to brightness value zero.
    * __Soft Hillshade__ -- Produces a hillshade layer with shade values as floating point numbers from zero to one. 
      This output includes the entire range of output from Lambert (-1 to 1), but re-scales it to a 'normalized'
      range from zero to one.  In this output, values between 0 and 0.5 represent 'hidden' or 'self-shaded' areas
      where the surface normal vector is greater than 90 degrees from the light vector. Shadows are also modeled.  
      Rather than casting shadows to full-black, a pixel found to be in shadow is multiplied by 0.5.  This will darken
      the pixel substantially, but permits variation among shadowed areas, preserving detail.  Recommend stretching 
      this output on a minimum-maximum color ramp where:
        * Hillshade values from 0 to 0.5 ramp from 100% black to 80% black
        * Hillshade values from 0.5 to 1.0 ramp from 80% black to 0% black (100% white)
      
        This will allow for subtle variation int the darkest parts of the scene, while providing ample range between
        80% black and 0% black for the majority of the scene's surface. 
        A style file is included in the repository which defines an '80-20' and a '70-30' color ramp for soft shaded
        scenes using this tool.  Ensure your stretch type is 'Min-Max' with the minimum set to zero and the maximum 
        set to one. 
    * __Sky Model__ -- Massively multi-directional hillshade; several hundred light sources supported. 
      This implementation is based on a hard hillshade using output cast as 8-bit integers. 
* _<u>NormalMaps</u>_ ToolSet
    * __NLCD_Normal_Mapper__ -- Configures and applies bump-map tiles according to a land use raster mask. The 
      use-case was built around land cover, which is why *NLCD* terminology is used. The bump-mapper expects 
      the NLCD masking raster to have same cell size, extent, dimensions, and projection as the terrain raster.
    * __Prep_NLCD_Normalmap_Mask__ -- This tool helps to clip, resample, and re-project NLCD or similar raster data to
      match the surface normal dataset.
      

Options for all Surface Normal Tools:
* Near the top of the `SurfaceOrientationToolbox.pyt` file is a declaration for global variable `DEBUG`. If this is
  set to <u>True</u>, debug information is accumulated in a log file in the same folder/directory
  as the toolbox source file.
  
### Normal Maps
The bump-mapping tool searches the `NormalMaps` folder for the tiles it offers for you to use in texturing 
a scene. You can add normal maps to this foldor *prior* to running the tool to give yourself more options. 
You can create your own bump-map tiles, or download most any normal map from on-line archives. 
The libraries at https://polyhaven.com/textures and https://ambientcg.com/ are popular repositories of 
public-domain textures. The Normal-Mapper tool recognizes tiff files only (and must be named with .tiff extension).

### Sky Config Files
The configuration files in the `SkyConfigFiles` folder provide the options for sky configurations in the
Sky Model tool.  If you want to create your own sky configurations, you can add them as a text file in this
folder to have them included as an option for the tool. 

### Jupyter Notebooks
The `SketchNotebooks` folder includes a few Jupyter notebooks which capture some of the algorithms
used in the toolbox.  Some of these use the Crater Lake DEM (as GeoTIFF) dataset.  You'll need to 
download your own copy of that dataset from https://shadedrelief.com/SampleElevationModels/ and adjust
the path to the TIFF in the notebook. 

### Support Files
* __`LambertColorRamps.stylx`__ -- ArcGIS Pro style file for the color ramps used for soft hillshades. 
* __`NLCD.clr`__ -- Color definitions for National Land Cover Database rasters.
* __`LCCS.clr`__ -- Color definitions for Copernicus land cover rasters.
