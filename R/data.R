#' A data frame of vegetation parameters for CORINE land cover types
#' @format a data frame with the following columns:
#' \describe{
#'  \item{Code}{Numeric code used in CORINE landcover data}
#'  \item{Descriptor}{Description of vegetation type}
#'  \item{x}{Ratio of vertical to horizontal projections of leaf foliage}
#'  \item{gsmax}{Maximum stomatal conductances (mol / m^2 / s)}
#'  \item{leafd}{Leaf diameter (m)}
#' }
"corinetable"

#' A template data frame for storing user credentials
#'
#' A data frame to be used as a template for storing user credentials
#'
#' @format a data frame with the following elements:
#' \describe{
#'  \item{Site}{Repository for which user credentials are being provided}
#'  \item{username}{username}
#'  \item{password}{password}
#'  \item{Datasets}{Descriptor of datasets accessed from that repository}
#'  \item{Weblink}{Weblink for registration}
#' }
"credentials"

#' A data frame of vegetation parameters for ESA land cover types
#' @format a data frame with the following columns:
#' \describe{
#'  \item{Code}{Numeric code used in ESA landcover data}
#'  \item{Descriptor}{Description of vegetation type}
#'  \item{x}{Ratio of vertical to horizontal projections of leaf foliage}
#'  \item{gsmax}{Maximum stomatal conductances (mol / m^2 / s)}
#'  \item{leafd}{Leaf diameter (m)}
#' }
"esatable"

#' A dataset of leaf area index values
#'
#' A wrapped SpatRaster of leaf area index values for the Lizard Peninsula, Cornwall,
#' UK, for February 2023 for the area bounded by 160000, 165000, 25000, 30000
#' (OSGB36 / British National Grid, EPSG:27700)
#' @format A PackedSpatRaster object with 500 rows and 500 columns
"laifeb"

#' NCEP-NOAA land-sea mask
#'
#' A wrapped SpatRaster indicating whether NCEP-NOAA grid cells are classed as
#' land (1) or sea (0) for the area bounded by -180.9375, 179.0625, -89.49406
#' (xmin, xmax, ymin, ymax) lon/lat WGS 84
#' @format A PackedSpatRaster object with 94 rows and 192 columns
#' "landsea"

#' NCEP-NOAA land-sea fraction
#'
#' A wrapped SpatRaster indicating what fraction of NCEP-NOAA grid cells
#' comprise land for the area bounded by -180.9375, 179.0625, -89.49406
#' (xmin, xmax, ymin, ymax) lon/lat WGS 84
#' @format A PackedSpatRaster object with 94 rows and 192 columns
#' "landseafrac"

#' A dataset of monthly leaf area index values
#'
#' A wrapped multi-layer SpatRaster of monthly leaf area index values for the Lizard Peninsula, Cornwall,
#' UK, for Jan-Dev 2023 for the area bounded by 160000, 162000, 28000, 30000
#' (OSGB36 / British National Grid, EPSG:27700)
#' @format A PackedSpatRaster object with 200 rows, 200 columns
"monthlylai"

#' A dataset of soilparameters parameters.
#'
#' An object of class groundparams
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{Soil.type}{description of soil type}
#'   \item{Smax}{Volumetric water content at saturation (m^3 / m^3)}
#'   \item{Smin}{Residual water content (m^3 / m^3)}
#'   \item{alpha}{Shape parameter of the van Genuchten model (cm^-1)}
#'   \item{n}{Pore size distribution parameter (dimensionless, > 1)}
#'   \item{Ksat}{Saturated hydraulic conductivity (cm / day)}
#'   \item{Vq}{Volumetric quartz fraction of soil}
#'   \item{Vm}{Volumetric mineral fraction of soil}
#'   \item{Vo}{Volumetric organic fraction of soil}
#'   \item{Mc}{Mass fraction of clay}
#'   \item{rho}{Soil bulk density (Mg / m^3)}
#'   \item{b}{Shape parameter for Campbell model (dimensionless, > 1)}
#'   \item{psi_e}{Matric potential (J / m^3)}
#'   \item{mult}{soil moisture model radiation coefficient}
#'   \item{rmu}{soil moisture model rainfall coefficient}
#'   \item{a}{soil moisture model deeper layer multiplier coefficient}
#'   \item{pwr}{soil moisture model deeper layer power coefficient}
#' }
#' @source: \url{https://onlinelibrary.wiley.com/doi/full/10.1002/ird.1751}
"soiltable"

#' A dataset of coastal exposure values for the UK
#'
#' A one km grid resolution wrapped multi-layer SpatRaster of coastal exposure values in 8 compass directions
#' (clockwise form north) for the area bounded by 0, 700000, 0, 1250000 (OSGB36 / British National Grid, EPSG:27700)
#' @format A PackedSpatRaster object with 1250 rows, 700 columns and 8 layers.
#' "ukcoastexposure"
