#' @title Download Landcover data
#'
#' @description Downloads either 10m grid resolution ESA or 100m grid resolution
#' CORINE landcover data from Google Earth Engine. The latter contains more habitat
#' classes.
#' @param r a SpatRaster with a coordinate reference system defined. Used
#' to the location for which to download data.
#' @param type one of `ESA` or `CORINE`
#' @param year year for which land cover data are required. Ignored if `type` set
#' to `ESA` as 2023 data will be downloaded. Options are 1986 to 2018 for `CORINE`.
#' @param GoogleDrivefolder name of folder on your Google drive to which data will
#' be downloaded.
#' @param pathtopython directory in which Python is saved (see
#' vignette instructions for installing and setting up Python)
#' @param projectname the name of your Google Earth Engine project (see
#' vignette instructions.
#' @param silent optional logical indicating whether or not to inform user of download
#' progress.
#' @details See vignette for example of how to use ths function.
#' @import terra
#' @import reticulate
#' @import rgee
#' @export
#' @rdname lcover_download
lcover_download<-function(r, type = "ESA", year = 2018, GoogleDrivefolder, pathtopython, projectname = NA, silent = FALSE) {
  typecheck <- type %in% c("ESA", "CORINE")
  if (typecheck == FALSE) stop("type must be one of ESA or CORINE")
  if (type == "CORINE") {
    yearcheck <- year %in% c(1986:2018)
    if (yearcheck == FALSE) stop("CORINE land cover data only available on Google Earth Engine for 1986 to 2018")
  }
  use_python(paste0(pathtopython,"python.exe"), required = TRUE)
  if (is.na(projectname) == FALSE) ee$Initialize(project=projectname)
  # Get bounding box
  e<-ext(r)
  r2<-rast(e)
  crs(r2)<-crs(r)
  r2<-project(r2,"EPSG:4326")
  e<-ext(r2)
  # work out crs
  proj_string <- crs(r, describe=T)
  epsg_code <- paste0("EPSG:",proj_string$code)
  # get aoi_coordinates
  aoi <- ee$Geometry$Rectangle(c(e$xmin, e$ymin, e$xmax, e$ymax))
  aoi_bounds <- aoi$bounds()$getInfo()
  aoi_coordinates <- aoi_bounds$coordinates[[1]]
  if (type == "ESA") {
    collection <- ee$ImageCollection('ESA/WorldCover/v100')
    lcover <- collection$first()
    task <- ee$batch$Export$image$toDrive(
      image = lcover,  # The image to export
      description = 'ESA_WorldCover_Export',  # Description for the export task
      folder = GoogleDrivefolder,  # Folder in Google Drive where the file will be saved
      fileNamePrefix = 'ESA_WorldCover_2023',  # Prefix for the filename
      region =   aoi_coordinates,  # Define the export region
      scale = 10,  # Scale in meters (resolution)
      crs = epsg_code  # Coordinate reference system
    )
  } else {
    # Define the CORINE Landcover dataset
    dataset <- ee$Image(paste0('COPERNICUS/CORINE/V20/100m/',year))
    # Select the 'landcover' band
    landCover <- dataset$select('landcover')
    task <- ee$batch$Export$image$toDrive(
      image = landCover,  # The image to export
      description = 'CORINE_LandCover_Export',  # Description for the export task
      folder = GoogleDrivefolder,  # Folder in Google Drive where the file will be saved
      fileNamePrefix = paste0('CORINE_LandCover_',year),  # Prefix for the filename
      region = aoi_coordinates,  # Define the export region
      scale = 100,  # Scale in meters (resolution) - for CORINE dataset
      crs = epsg_code   # Coordinate reference system
    )
  }
  task$start()
  if (silent == FALSE) .monitor_task(task$id)
}
#' @title Download Landcover data
#'
#' @description Downloads either 10m grid resolution ESA or 100m grid resolution
#' CORINE landcover data from Google Earth Engine. The latter contains more habitat
#' classes.
#' @param r a SpatRaster with a coordinate reference system defined. Used
#' to the location for which to download data.
#' @param GoogleDrivefolder name of folder on your Google drive to which data will
#' be downloaded.
#' @param pathtopython directory in which Python is saved (see
#' vignette instructions for installing and setting up Python)
#' @param projectname the name of your Google Earth Engine project (see
#' vignette instructions).
#' @param silent optional logical indicating whether or not to inform user of download
#' progress.
#' @details See vignette for example of how to use ths function.
#' @import terra
#' @import reticulate
#' @import rgee
#' @export
#' @rdname vegheight_download
vegheight_download<-function(r, GoogleDrivefolder, pathtopython, projectname = NA, silent = FALSE) {
  use_python(paste0(pathtopython,"python.exe"), required = TRUE)
  if (is.na(projectname) == FALSE) ee$Initialize(project=projectname)
  # Get bounding box
  e<-ext(r)
  r2<-rast(e)
  crs(r2)<-crs(r)
  r2<-project(r2,"EPSG:4326")
  e<-ext(r2)
  # work out crs
  proj_string <- crs(r, describe=T)
  epsg_code <- paste0("EPSG:",proj_string$code)
  # get aoi_coordinates
  aoi <- ee$Geometry$Rectangle(c(e$xmin, e$ymin, e$xmax, e$ymax))
  aoi_bounds <- aoi$bounds()$getInfo()
  aoi_coordinates <- aoi_bounds$coordinates[[1]]
  canopy_height <- ee$Image('users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1')
  task <- ee$batch$Export$image$toDrive(
    image = canopy_height,  # The image to export
    description = 'canopy_height_export',  # Description for the export task
    folder = GoogleDrivefolder,  # Folder in Google Drive where the file will be saved
    fileNamePrefix = 'canopy_height_2020',  # Prefix for the filename
    region =   aoi_coordinates,  # Define the export region
    scale = 10,  # Scale in meters (resolution)
    crs = epsg_code  # Coordinate reference system
  )
  task$start()
  if (silent == FALSE) .monitor_task(task$id)
}
#' @title Download Leaf Area Index (LAI) data
#'
#' @description Downloads either 10m grid resolution LAI data from WEkEO or
#' 500m grid resolution LAI data from LPDAAC_ECS.
#' @param r a SpatRaster with a coordinate reference system defined. Used
#' to the location for which to download data.
#' @param tme POSIXlt object covering time period for which data are required
#' (see details).
#' @param reso the spatial resolution of data to download. Must be either
#' 10 or 500 (see details).
#' @param pathout directory to which data will be downloaded.
#' @param credentials a data.frame of user credentials (see vignette for
#' setting this up)
#' @param pathtopython directory in which Python is saved (see
#' vignette instructions for installing and setting up Python)
#' @seealso [lai_mosaic()], [lai_fill()], [lai_downscale()],
#' [lai_seasonadj()]
#' @details if `reso = 10` the R package `reticulate` is used to interface
#' with Python and download 10m resolution LAI data from the Copernicus WEkEO
#' site. Data are only available from 2019 to present and only for parts of
#' Europe (see https://www.wekeo.eu/data?view=dataset&dataset=EO%3AEEA%3ADAT%3ACLMS_HRVPP_VI
#' for details of coverage).If `reso = 500` 500 m grid resolution LAI data derived
#' from MODIS imagery are downloaded from LPDAAC_ECS using the `luna` package.
#' Data are available globally from 18th February 2018 onward. Multiple
#' files are downloaded and stored in the directory given by `pathout`. Use
#' [lai_mosaic()] to tile these together and [lai_fill()] to fill data gaps
#' caused by cloud cover. See vignette for examples.
#' @import luna
#' @import terra
#' @import reticulate
#' @export
#' @rdname lai_download
lai_download<-function(r, tme, reso = 10, pathout, credentials, pathtopython = "C:/Python/") {
  isres <- reso %in% c(10, 500)
  if (isres == FALSE) stop("reso must be one of 10 or 500")
  year<-tme$year[1]+1900
  if (reso == 10) {
    cyear<-as.POSIXlt(Sys.time())$year+1900
    isyr <- year %in% c(2019:cyear)
    if (isyr == FALSE) stop("Data not available for specified time period")
    # Configure user's credentials if .hdarc not existing
    use_python(paste0(pathtopython,"python.exe"))
    # Import HDA Python module
    hda <- import("hda")
    username <- credentials$username[2]
    password <- credentials$password[2]
    conf <- hda$Configuration(user =  username, password = password)
    hda_client <- hda$Client(config = conf)
    # create bounding box
    bb<-.bbox(r)
    # divide tme into 5 day time slices
    tmes<-seq(tme[1],tme[length(tme)],by=5*3600*24)
    ttxt<-paste0(substr(tmes,1,10),"T00:00:00.000Z")
    n<-length(ttxt)
    st<-ttxt[1:(n-1)]
    ed<-ttxt[2:n]
    # Do stuff for 10 m res
    files<-0
    for (i in 1:length(st)) {
      # Set up query for LAI data
      query <- list(dataset_id = "EO:EEA:DAT:CLMS_HRVPP_VI",
                    productType = "LAI",
                    resolution = "10",
                    start = st[i],
                    end = ed[i],
                    bbox = bb)
      matches <- hda_client$search(query)
      if (length(matches) > 0) {
        files<-files+length(matches)
        matches$download(pathout)
        # Set up query for Quality control data
        query <- list(dataset_id = "EO:EEA:DAT:CLMS_HRVPP_VI",
                      productType = "QFLAG2",
                      resolution = "10",
                      start = st[i],
                      end = ed[i],
                      bbox = bb)
        matches <- hda_client$search(query)
        matches$download(pathout)
      }
    }
    if(files == 0) stop("No data for specified location or time period")
  } else {
    tsed<-tme[length(tme)]
    tmod<-as.POSIXlt(0,origin="2000-02-18",tz="UTC")
    if (tsed < tmod) stop("No data available prior to 2000-02-18")
    # Download MODIS DATA
    e<-ext(r)
    r2<-rast(e)
    crs(r2)<-crs(r)
    r2<-project(r2, "EPSG:4326")
    e2<-ext(r2)
    st<-substr(as.character(tme[1]),1,10)
    ed<-substr(as.character(tme[length(tme)]),1,10)
    mf<-luna::getNASA("MOD15A2H", st, ed, aoi = e2, version = "061", download=FALSE)
    username <- credentials$username[1]
    password <- credentials$password[1]
    if (length(mf) > 0) {
      luna::getNASA("MOD15A2H", st, ed, aoi = e2, version = "061", download=TRUE,
                    path=pathout,username=username,password=password,server="LPDAAC_ECS")
    } else {
      stop("No data for specified location or time period")
    }
  }
}
#' @title Moasic Leaf Area Index (LAI) data returned
#'
#' @description This function moasics Leaf Area Index (LAI) data downloaded using
#' `lai_download`.
#'
#' @param r a SpatRaster with a coordinate reference system defined. Used
#' for masking and reprojecting the downloaded data.
#' @param pathin filepath of downloaded data. Same as `pathout` in function `lai_download`.
#' @param reso the spatial resolution of data download. Must be either
#' 10 or 500 (see details under `lai_download`).
#' @param msk option logical indicating whether data should be masked using `r`.
#' @seealso [lai_download()], [lai_fill()]
#' @import terra
#' @export
#' @rdname lai_mosaic
lai_mosaic <- function(r, pathin, reso = 10, msk = TRUE) {
  isres <- reso %in% c(10, 500)
  if (isres == FALSE) stop("reso must be one of 10 or 500")
  # Create list of files
  lst<-list.files(pathin)
  if (reso == 10) {
    wlai<-rep(0,length(lst))
    # Select only those that are LAI
    for (i in 1:length(lst)) {
      n<-nchar(lst[[i]])
      ed<-n-4
      st<-n-6
      ss<-substr(lst[[i]],st,ed)
      if (ss=="LAI") wlai[i]<-1
    }
    s<-which(wlai==1)
    lfile<-lst[s]
    # Create list of qflag files
    qfile<-paste0(substr(lfile,1,40),"QFLAG2.tif")
    # Loop through each and:
    ro<-list()
    tile<-0
    n<-length(lfile)
    pb<-utils::txtProgressBar(min = 0, max = n*1.5, style = 3)
    for (i in 1:n) {
      utils::setTxtProgressBar(pb, i)
      # (a) check for spatial overlap
      fi<-paste0(pathin,lfile[i])
      ri<-rast(fi)
      e<-ext(ri)
      re<-rast(e)
      crs(re)<-crs(ri)
      if (crs(r) != crs(re)) re<-project(re,crs(r))
      oip<-intersect(ext(r),ext(re))
      if (is.null(oip) == FALSE) {
        xx<-suppressWarnings(.chcktif(fi))
        if (is.na(xx) == FALSE) {
          fiq<-paste0(pathin,qfile[i])
          if (file.exists(fiq)) {
            qu<-rast(fiq)
            ri<-.cleanlai(ri,qu)
            if (crs(ri) != crs(r)) ri<-project(ri,r)
            ri<-resample(ri,r,method="near")
            if (msk) ri<-mask(ri,r)
            v<-as.vector(r)
            v<-v[is.na(v) == FALSE]
            if (length(v) > 0) {
              tile<-tile+1
              ro[[tile]]<-ri
            } # end NA check
          } # end quality file exists check
        } # end file corrupted check
      } # end overlap check
    } # end loop
    # mosaic
    rma<-ro[[1]]
    nn<-length(ro)
    if (length(ro) > 1) {
      for (i in 2:length(ro)) {
        utils::setTxtProgressBar(pb, ((i*0.5*n)/nn)+n)
        rr<-ro[[i]]
        rma<-mosaic(rma,rr)
      } # end for
    } # end if
  } # end reso
  else {
    ro<-list()
    tile<-0
    n<-length(lst)
    pb<-utils::txtProgressBar(min = 0, max = n, style = 3)
    for (i in 1:n) {
      utils::setTxtProgressBar(pb, i)
      fi<-paste0(pathin,lst[i])
      ri<-rast(fi)[[2]]
      mskr<-rast(fi)[[4]]
      ri<-mask(ri,mskr)
      if (crs(ri)!=crs(r)) {
        e<-ext(r)
        r2<-rast(e)
        crs(r2)<-crs(r)
        r2<-project(r2,crs(ri))
        oip<-intersect(ext(r2),ext(ri))
        if (is.null(oip) == FALSE) {
          tile<-tile+1
          ri<-crop(ri,ext(r2))
          ri<-project(ri,crs(r))
          rr<-rast(ext(ri))
          res(rr)<-500
          crs(rr)<-crs(ri)
          ro[[tile]]<-resample(ri,rr)
        }
      } else {
        oip<-intersect(ext(r),ext(ri))
        if (is.null(oip) == FALSE) {
          tile<-tile+1
          ri<-crop(ri,ext(r))
          rr<-rast(ext(ri))
          res(rr)<-500
          crs(rr)<-crs(ri)
          ro[[tile]]<-resample(ri,rr)
        }
      }
    }
    rma<-ro[[1]]
    if (length(ro) > 1) {
      for (i in 2:length(ro)) {
        rr<-ro[[i]]
        rma<-mosaic(rma,rr)
      } # end for
    }
    if (msk) {
      rmsk<-resample(r,rma)
      rma<-mask(rma,rmsk)
    }
  }
  return(rma)
}
#' @title Cloud cover fill lai data using landcover data
#'
#' @description This function cover fill Leaf Area Index (LAI) using land cover data
#'
#' @param lai a SpatRaster of Leaf Area Index values with some NAs
#' @param lcover a SpatRaster of land cover matching the resolution, extent and
#' coordinate reference system of `lai`
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimdata, .registration = TRUE
#' @import terra
#' @export
#' @rdname lai_fill
#' @return a cloud filled SpatRaster of Leaf Area Index values.
lai_fill<-function(lai, lcover) {
  lai<-mask(lai,lcover)
  n<-dim(lai)[3]
  if (n == 1) {
    if (crs(lai) != crs(lai)) stop("projections of lai and lcover must match")
    if (dim(lai)[1] != dim(lcover)[1]) stop("grid cells of lai and lcover must match")
    if (dim(lai)[2] != dim(lcover)[2]) stop("grid cells of lai and lcover must match")
    lai<-fillpai(.is(lai),.is(lcover))
    lai<-.rast(lai,lcover)
  } else {
    ## create lpai matrix
    mpaif<-function(lai,lcover,u,i) {
      mpai<-meanpai(.is(lai[[i]]),.is(lcover))
      if (length(mpai) > 0) {
        ss<-seq(min(u,na.rm=TRUE),max(u,na.rm=TRUE))
        matches <- which(ss %in% u)
        ss<-ss*NA
        ss[matches]<-mpai
      } else ss<-NA
      ss
    }
    pb<-utils::txtProgressBar(min = 0, max = n+3, style = 3)
    u<-uniquecpp(.is(lcover))
    u<-u[is.na(u)==FALSE]
    utils::setTxtProgressBar(pb, 1)
    ss<-mpaif(lai,lcover,u,1)
    lpai<-matrix(ss,ncol=n,nrow=length(ss))
    for (i in 2:n) {
      utils::setTxtProgressBar(pb, i)
      ss<-mpaif(lai,lcover,u,i)
      lpai[,i]<-ss
    }
    for (i in 1:dim(lpai)[1]) {
      me<-mean(lpai[i,],na.rm=TRUE)
      if (is.na(me) == FALSE) lpai[i,]<-na_approx(lpai[i,])
    }
    lai<-napproxCpp(.is(lai), .is(lcover), lpai)
    utils::setTxtProgressBar(pb, n+3)
    lai<-.rast(lai,lcover)
  }
  return(lai)
}
#' @title Downscales LAI data
#'
#' @description This function uses high-resolution land cover and/or vegetation
#' height data to downscale coarser resolution LAI.
#'
#' @param lai a coarse-resolution SpatRaster of Leaf Area Index values
#' @param lcover a high-resolution SpatRaster of land cover matching the
#' coordinate reference system of `lai`
#' @param veghgt optionally, a high-resolution SpatRaster of vegetation
#' heights matching the resolution, extent  and coordinate reference system of `lcover`
#' @param aszero optionally a vector of values indicating, for which land cover types,
#' `lai` values should automatically be set at zero
#' @import terra
#' @export
#' @rdname lai_downscale
#' @return a downscaled SpatRaster of Leaf Area Index values.
#' @details if `veghgt` is not provided only land cover data are used to perform
#' the downscaling.
lai_downscale<-function(lai, lcover, veghgt = NA, aszero = NA) {
  if (class(veghgt) != "logical") {
    veghgt[is.na(veghgt)]<-0
    veghgtc<-resample(veghgt,lai)
    y<-as.vector(veghgtc)
    x<-as.vector(lai)
    s<-which(y>0 & x>0)
    y<-y[s]
    x<-x[s]
    m1<-summary(lm(log(y)~log(x)))
    epai<-m1$coef[1,1]+m1$coef[2,1]*log(.is(veghgt))
    epai<-.rast(exp(epai),veghgt)
    epai[veghgt==0]<-0
    epaic<-resample(resample(epai,lai),veghgt)
    epai<-suppressWarnings(epai/epaic)
    epai[epai>5]<-5
    epai[is.na(epai)]<-0
  }  else {
    epai<-lcover*0+1
  }
  if (class(lcover) != "logical") {
    # Produce table of habitat types
    lcover<-lcover+1 # add one to ensure there are no zeros
    habs<-unique(as.vector(lcover))
    s<-which(is.na(habs))
    if (length(s) > 0) {
      nas<-max(habs,na.rm=T)+1
      habs[.is.na(habs)]<-nas
      habs<-unique(as.vector(lcover))
    }
    mlai<-0
    for (i in 1:length(habs)) {
      # Calculate fraction of each habitat type in each lai grid cell
      htype<-lcover
      htype[htype!=habs[i]]<-0
      htype[htype==habs[i]]<-1
      habfrac<-resample(htype,lai)
      # multiply LAI by fraction and sum
      mulai<-habfrac*lai
      hsum<-sum(as.vector(mulai),na.rm=TRUE)
      # Sum the habitat type fractions across cells
      fsum<-sum(as.vector(habfrac),na.rm=TRUE)
      # Calculate the mean LAI for each habitat type
      mlai[i]<-hsum/fsum
    }
    epai2<-.is(lcover)*0
    for (i in 1:length(habs)) {
      s<-which(.is(lcover)==habs[i])
      epai2[s]<-mlai[i]
    }
    epai2<-.rast(epai2,lcover)
    epai<-suppressWarnings(epai2*epai)
    epai[is.na(epai)]<-0
  }
  # resacle so coarse grid cell average is one
  af<-res(lai)[1]/res(lcover)[1]
  cpai<-aggregate(epai,af)
  cpai<-resample(cpai,lcover)
  mupai<-epai/cpai
  mupai[is.na(mupai)]<-0
  flai<-resample(lai,lcover)
  flai[is.na(flai)]<-0
  flai<-.is(flai*mupai)
  # Cap at twice epai
  mxpai<-.is(2*epai)
  s<-which(flai>mxpai)
  flai[s]<-mxpai[s]
  if (class(aszero) != "logical") {
    for (i in 1:length(aszero)) {
      s<-which(.is(lcover) == aszero[i]+1)
      flai[s]<-0
    }
  }
  if (class(veghgt) != "logical") {
    r<-veghgt
  } else r<-lcover
  flai<-.rast(flai,r)
  if (class(veghgt) != "logical") flai[veghgt==0]<-0
  mx<-max(as.vector(lai),na.rm=TRUE)*2
  flai[flai>mx]<-mx
  flai
}
#' @title Cloud fills a multi-layer SpatRaster of LAI values by applying season effects
#'
#' @description This function takes as an input a multi-layer SpatRaster of e.g. monthly
#' LAI values for a year. Missing values are interpolated temporally.
#'
#' @param lai a multi-layer SpatRaster of Leaf Area Index values
#' with some missing values
#' @returns a multi-layer SpatRaster of Leaf Area Index values
#' without missing values
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimdata, .registration = TRUE
#' @import terra
#' @export
lai_seasonadj<-function(lai) {
  nms<-names(lai)
  laia<-as.array(lai)
  seffect<-seasoneffect(laia)
  # Check whether any NAs in season effect
  m<-mean(seffect)
  if (is.na(m)) {
    if (is.na(seffect[1])) {
      seffect<-.circular_approx(seffect)
    } else {
      seffect<-na_approx(seffect)
    }
  }
  laia<-seasonadjCpp(laia,seffect)
  lai<-.rast(laia,lai[[1]])
  names(lai)<-nms
  return(lai)
}
#' @title Derive leaf area index from NDVI
#'
#' @description This function derives leaf area index data from an estimate of the
#' normalized difference vegetation index (NDVI).
#'
#' @param RGB a red-green-blue image supplied as a SpatRaster (see details).
#' @param CIR a colour-infrared image supplied as a SpatRaster (see details).
#' @param modisLAI a coarser resolution SpatRaster of MODIS-derived lai data
#' (see [lai_download()]).
#' @param maxlai optional maximum leaf area index value expected within study (see details)
#' @returns a SpatRaster of Leaf Area Index values
#' @details `RGB` and `CIR` are normally 3-band multi-layer SpatRasters derived from
#' aerial photographs. It is assumed that the first band in the `RGB` image corresponds
#' to red and the first band in the `CIR` corresponds to near-infrared (NIR). NDVI
#' is then calculated as NDVI = (NIR - RED) / (NIR + RED). Since only the first
#' bands of `RGB` and `CIR` are used for this calaculation, `RGB` and `CIR` can
#' instead be single layer SpatRasters of red and near-infrared reflectance respectively.
#' Since aerial photographs are normally brightness and contrast adjusted during processing,
#' the derived NDVI values are not true NDVI. Thus the `modisLAI` SpatRaster is used to
#' rescale derived lai. `maxlai` is also used for rescaling, but if not supplied, is
#' estimated from the maximum value in `modisLAI` adjusted upwards depending on the resolution
#' of `RGB` and `CIR` relative to `modisLAI`. `RGB` and `CRS` must have matching
#' coordinate reference systems, extents and resolutions, whereas `modisLAI` is typically
#' of coarser resolution and reprojected if necessary.
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimdata, .registration = TRUE
#' @import terra
#' @export
lai_fromndvi <- function(RGB, CIR, modisLAI, maxlai = NA) {
  # calculate ndvi
  ndvi <- (CIR[[1]] - RGB[[1]]) / (CIR[[1]] + RGB[[1]])
  # find maximum pai
  if (class(maxlai) == "logical") {
    af <- res(modisLAI)[1] / res(RGB)[1]
    mu <- 1 / (0.4468498 + 0.0064626 * af)
    mu[mu < 1] <- 1
    maxlai <- max(as.vector(modisLAI), na.rm = TRUE) * mu
  }
  # find max ndvi for rescaling
  S1 <- exp(-0.956782804 * maxlai)
  D1 <- -0.682265131 / S1 + 0.00796106 * S1
  p1 <- (0.061666667 / (D1 * S1)) * -0.357497089
  p2 <- (-0.061666667 * S1 / D1) * 1.556068518
  albr <- p1 + p2
  mxndvi <- (0.212278228 - albr) / (0.212278228 + albr)
  # rescale ndvi
  mu <- (mxndvi - 0.096256966) / (max(as.vector(ndvi), na.rm = TRUE) - min(as.vector(ndvi), na.rm = TRUE))
  ndvi <- ndvi * mu
  ndvi <- ndvi - min(as.vector(ndvi), na.rm = TRUE) + 0.096256966
  # find pai
  pai <- .rast(find_pai(as.matrix(ndvi, wide = TRUE)), ndvi)
  # adjust using modis layer
  if (crs(modisLAI) != crs(ndvi)) modisLAI <- project(modisLAI, crs(ndvi))
  mu <- (modisLAI + 0.001) / (resample(pai, modisLAI) + 0.001)
  mu[is.na(modisLAI)] <- 1
  mu <- resample(mu, pai)
  mu[is.na(pai)] <- 1
  mu[is.na(mu)] <- 1
  pai <- pai * mu
  pai[pai > maxlai] <- maxlai
  pai <- mask(pai, ndvi)
  return(pai)
}
#' @title Calculate leaf inclination coefficient
#'
#' @description Estim,ates typical values of the ratio of vertical to horizontal
#' projections of leaf foliage from land cover data,
#'
#' @param landcover a SpatRaster of landcover class as returned by [lcover_download()].
#' @param lctype one of `ESA`, `CORINE`, or `Copernicus` indicating the source of `landcover`.
#' @returns a SpatRaster leaf inclination coefficient values
#' @import terra
#' @export
#' @rdname x_calc
x_calc <- function(landcover, lctype = "ESA") {
  u <- unique(as.vector(landcover))
  u <- u[is.na(u) == FALSE]
  x <- landcover * 0 + 0.01
  if (lctype == "ESA") {
    ltable <- esatable
  } else if (lctype == "CORINE") {
    ltable <- corinetable
  } else if (lctype == "Copernicus") {
    ltable <- gdlctable
  }
  for (i in 1:length(u)) {
    s <- which(ltable$Code == u[i])
    x[landcover == u[i]] <- ltable$x[s]
  }
  return(x)
}
#' @title Create vegetation inputs for grid model
#'
#' @description Creates an object of class vegparams used as an input to `microclimf`
#' or `microgrid`.
#'
#' @param landcover a SpatRaster of landcover class as returned by [lcover_download()].
#' @param vhgt a SpatRaster of vegetation height as returned by [vegheight_download()].
#' @param lai a SpatRaster of Leaf Area Index values.
#' @param refldata a list of ground and leaf reflectance values as returned by [reflectance_calc()].
#' @param lctype one of `ESA` or `CORINE` indicating the source of `landcover`. Alternatively a data.frame
#' matching vegetation parameters to landcover types following the format of the inbuilt datasets `esatable`
#' and `corinetable`.
#' @returns an object of class `vegparams`, namely a list of the following objects:
#' \describe{
#'   \item{pai}{a PackedSpatRaster of plant area index values}
#'   \item{hgt}{a PackedSpatRaster object of vegetation heights (m)}
#'   \item{x}{a PackedSpatRaster object of the ratio of vertical to horizontal projections of leaf foliage}
#'   \item{gsmax}{a PackedSpatRaster object maximum stomatal conductances (mol / m^2 / s)}
#'   \item{leafr}{a PackedSpatRaster object of leaf reflectance values}
#'   \item{clump}{a PackedSpatRaster of monthly values between 0 and 1 indicating the fraction of radiation passing through larger gaps in the canopy}
#'   \item{leafd}{a PackedSpatRaster object of leaf diameters (m)}
#'   \item{leaft}{a PackedSpatRaster object of leaf transmittance}
#' }
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimdata, .registration = TRUE
#' @import terra
#' @export
#' @rdname create_veggrid
create_veggrid <- function(landcover, vhgt, lai, refldata, lctype = "ESA") {
  # check geometries
  cp <- compareGeom(landcover, vhgt)
  if (cp) cp <- compareGeom(landcover, lai)
  if (cp) cp <- compareGeom(landcover, refldata$lref)
  msk <- landcover * vhgt * lai * refldata$lref
  landcover <- mask(landcover, msk)
  vhgt <- mask(vhgt, msk)
  lai <- mask(lai, msk)
  refldata$lref <- mask(refldata$lref, msk)
  # do lookups from land cover
  u <- unique(as.vector(landcover))
  u <- u[is.na(u) == FALSE]
  x <- gsmax <- leafd <- vhgt * 0 + 0.01
  if (lctype == "ESA") {
    ltable <- esatable
  } else if (lctype == "CORINE") {
    ltable <- corinetable
  } else if (lctype == "Copernicus") {
    ltable <- gdlctable
  } else {ltable <- lctype}
  for (i in 1:length(u)) {
    s <- which(ltable$Code == u[i])
    x[landcover == u[i]] <- ltable$x[s]
    gsmax[landcover == u[i]] <- ltable$gsmax[s]
    leafd[landcover == u[i]] <- ltable$leafd[s]
  }
  # estimate clumping factor
  nlyrs <- dim(lai)[3]
  clump <- lai * 0
  for (i in 1:nlyrs) {
    clump[[i]] <- .rast(calcclumpcpp(as.matrix(leafd, wide=TRUE),
                                     as.matrix(vhgt, wide=TRUE),
                                     as.matrix(lai[[i]], wide=TRUE)), x)
  }
  leafr <- refldata$lref
  leaft <- refldata$lref / 3
  # names
  names(lai) <- "Plant area index"
  names(vhgt) <- "Vegetation height"
  names(x) <- "Leaf inclination coefficient"
  names(gsmax) <- "Maximum stomatal conductance"
  names(leafr) <- "Leaf reflectance"
  names(clump) <- "Clumping coefficient"
  names(leafd) <- "Leaf diameter"
  names(leaft) <- "Leaf transmittance"
  vegp <- list(pai = wrap(lai),
               hgt = wrap(vhgt),
               x = wrap(x),
               gsmax = wrap(gsmax),
               leafr = wrap(refldata$lref),
               clump = wrap(clump),
               leafd = wrap(leafd),
               leaft = wrap(refldata$lref / 3))
  class(vegp) <- "vegparams"
  return(vegp)
}
#' @title Create vegetation inputs for point model
#'
#' @description Creates an object of class vegparams used as an input to `micropoint`.

#' @param landcover a SpatRaster of landcover class as returned by [lcover_download()].
#' @param vhgt a SpatRaster of vegetation height as returned by [vegheight_download()].
#' @param lai a SpatRaster of Leaf Area Index values.
#' @param refldata a list of ground and leaf reflectance values as returned by [reflectance_calc()].
#' @param lctype one of `ESA` or `CORINE` indicating the source of `landcover`. Alternatively a data.frame
#' matching vegetation parameters to landcover types following the format of the inbuilt datasets `esatable`
#' and `corinetable`.
#' @param lat optionally, latitude (decimal degrees) - see details.
#' @param long optionally, longitude (decimal degrees) - see details.
#' @details If `lat` and `long` are unspecified, point data are selelected for
#' the location and the centre of supplied `landcover`.
#' @returns an object of class `vegparams`, namely a list of the following objects:
#' \describe{
#'   \item{pai}{Plant area index value}
#'   \item{hgt}{Vegetation height (m)}
#'   \item{x}{Ratio of vertical to horizontal projection of leaf foliage}
#'   \item{gsmax}{Maximum stomatal conductances (mol / m^2 / s)}
#'   \item{leafr}{Leaf reflectance}
#'   \item{clump}{Value between 0 and 1 indicating the fraction of radiation passing through larger gaps in the canopy}
#'   \item{leafd}{Leaf diameter (m)}
#'   \item{leaft}{Leaf transmittance}
#' }
#' @import terra
#' @export
#' @rdname create_vegpoint
create_vegpoint <- function(landcover, vhgt, lai, refldata, lctype = "ESA", lat = NA, long = NA)  {
  if (class(lat) == "logical") {
    e<-ext(landcover)
    xy <- data.frame(x = (e$xmin + e$xmax) / 2,
                     y = (e$ymin + e$ymax) / 2)
  } else {
    # convert latlong to local crs
    xy <- data.frame(x = long, y = lat)
    xy <- sf::st_as_sf(xy, coords = c("x", "y"),
                       crs = 4326)
    xy <- sf::st_transform(xy, crs(landcover))
    xy <- data.frame(x = sf::st_coordinates(xy)[1], y = sf::st_coordinates(xy)[2])
  }
  reso <- res(landcover)[1]
  xmin <- floor((xy$x - reso) / reso) * reso
  xmax <- xmin + 2 * reso
  ymin <- floor((xy$y - reso) / reso) * reso
  ymax <- ymin + 2 * reso
  e <- ext(xmin, xmax, ymin, ymax)
  lcc <- crop(landcover, e)
  vhgtc <- crop(vhgt, e)
  laic <- crop(lai, e)
  refldata$lref <- crop(refldata$lref, e)
  vegp <- create_vegpgrid(lcc, vhgtc, laic, refldata, lctype)
  x <- rast(vegp$x)
  gsmax <- rast(vegp$gsmax)
  clump = rast(vegp$clump)
  leafd = rast(vegp$leafd)
  leaft = rast(vegp$leaft)
  vegpp <- list(h = extract(vhgtc, xy)[,2],
                pai = extract(laic, xy)[,2],
                x = extract(x, xy)[,2],
                clump = extract(clump, xy)[,2],
                lref = extract(refldata$lref, xy)[,2],
                ltra = extract(leaft, xy)[,2],
                leafd = extract(leafd, xy)[,2],
                em = 0.97,
                gsmax = extract(gsmax, xy)[,2],
                q50 = 100)
  class(vegpp) <- "vegparams"
  return(vegpp)
}
#' @title Create vegetation inputs for point model from create_veggrid output
#'
#' @description Creates an object of class vegparams used as an input to `micropoint`.

#' @param veggrid output of [create_veggrid()]. Can be unwrapped to save a step when processing many points
#' @param lat latitude
#' @param long longitude
#' @param llcrs the crs of lat,long coordinates
#' @param gediout optionally, single row of output from gedi_process for the point of interest.
#' Used to replace pai and canopy height outputs in veggrid
#' @returns an object of class `vegparams`, namely a list of the following objects:
#' \describe{
#'   \item{pai}{Plant area index value}
#'   \item{hgt}{Vegetation height (m)}
#'   \item{x}{Ratio of vertical to horizontal projection of leaf foliage}
#'   \item{gsmax}{Maximum stomatal conductances (mol / m^2 / s)}
#'   \item{leafr}{Leaf reflectance}
#'   \item{clump}{Value between 0 and 1 indicating the fraction of radiation passing through larger gaps in the canopy}
#'   \item{leafd}{Leaf diameter (m)}
#'   \item{leaft}{Leaf transmittance}
#' }
#' @import terra
#' @export
#' @rdname create_vegpoint
vegpfromgrid <- function(veggrid, lat, long, llcrs, pai, ch) {
  if(is(veggrid, "vegparams")) {
    veggrid = lapply(veggrid, unwrap)
  } 
  veggrid = rast(veggrid)
  pt = data.frame(lon=long, lat=lat)
  pt = vect(pt, crs = llcrs)
  pt = project(pt, veggrid[[1]])
  pt = terra::extract(veggrid, pt)
  vegpp <- list(h = pt$hgt,
                pai = pt$pai,
                x = pt$x,
                clump = pt$clump,
                lref = pt$leafr,
                ltra = pt$leaft,
                leafd = pt$leafd,
                em = 0.97,
                gsmax = pt$gsmax,
                q50 = 100)
  if(!is.null(gediout)) {
    vegpp$pai = gediout$pai
    vegpp$h = gediout$ch
  }
  class(vegpp) <- "vegparams"
  return(vegpp)
}

#' @title Process GEDI data for use in point-based microclimate models
#'
#' @description Extracts data from the L2B profile and calculates vertical PAI profile. The sum of the vertical PAI profile = total PAI

#' @param l2b character vector of hd5 files of GEDI L2B data (downloaded using gedi_download)
#' @param r spatRaster or spatVector to crop gedi data
#' @param sv optionally; provide a spatvector to further crop the gedi data
#' @param powerbeam default TRUE; if TRUE filters to only power beams; if FALSE includes power and coverage beams
#' @param yr optional; 4 digit year (YYYY) to filter GEDI data
#' @param mth optional; 2 digit month (MM) to filter GEDI data
#' @returns a dataframe of vegetation information with the following values
#' \describe{
#'   \item{pai}{Total plant area index value}
#'   \item{pai_z_}{plant area index between z1 to z2 m above ground}
#' }
#' @import rGEDI
#' @export
gedi_process<-function(l2b, r, sv, powerbeam=TRUE, yr=NULL, mth=NULL) {
  gedi = list()
  if(crs(r, proj=T)!= "+proj=longlat +datum=WGS84 +no_defs") r = project(r, "epsg:4326")
  for(i in 1:length(l2b)) {
    l2b_i = tryCatch(rGEDI::readLevel2B(l2b[i]), 
                     error = function(e){
                       message("GEDI file corrupt: redownloading")
                       return(NA)
                       })
    
    # redownload if the file is corrupt
    if(!is(l2b_i, "gedi.level2b")) {
      outdir = paste0(dirname(l2b[i]), "/")
      bname = basename(l2b[i])
      fpath = paste0("https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/GEDI02_B.002/",
                     gsub(".h5", "", bname), "/", bname)
      gediDownload(fpath, 
                   outdir = outdir,
                   overwrite = T)
      l2b_i = tryCatch(rGEDI::readLevel2B(l2b[i]), 
                       warning = function(e){error("GEDI file corrupt: redownloading")})
    }
    
    l2b_i = .getL2Bprofile(l2b_i)
    
    e = ext(r)
    e = project(e, from=crs(r), to="epsg:4326")
    
    # filter out poor quality shots and crop to extent in r
    l2b_i = l2b_i %>% 
      filter(l2b_quality_flag > 0) %>% 
      filter(lon_lowestmode >= e[1] & lon_lowestmode <= e[2] & lat_lowestmode >= e[3] & lat_lowestmode <= e[4]) %>% 
      # see l2b data dictionary for calculation of shot time
      mutate(shot_time = as.POSIXlt("2018-01-01") + delta_time,
             # convert rh100 from cm to m
             rh100 = rh100/100)
    
    if(!is.null(sv)) {
      l2b_i = l2b_i %>% 
        vect(geom = c("lon_lowestmode", "lat_lowestmode"), crs = "epsg:4326", keepgeom=T) %>% 
        crop(gedicrop) %>% 
        as.data.table()
    }
    
    
    # filter to power beams, year, and month if requested
    if(powerbeam) {
      l2b_i = l2b_i %>% 
        filter(beam == "BEAM0101" | beam == "BEAM0110" | beam == "BEAM1000" | beam == "BEAM1011")
    }
    if(!is.null(yr)) {
      l2b_i[, yr := year(shot_time)]
      l2b_i = l2b_i[yr==yr,]
    }
    if(!is.null(mth)) {
      l2b_i[, mth := month(shot_time)]
      l2b_i = l2b_i[mth==mth,]
    }
    
   
    # calculate vertical profile for use in run_micropoint. sum of the profile = total PAI (run check before running micropoint)
    # to prevent R crashing
    if(nrow(l2b_i)>0) {
      pai_cols <- paste0("pai_z", seq(0, 145, 5), "_", seq(5, 150, 5), "m")
      paiz = l2b_i[,.SD, .SDcols = pai_cols]
      paiz = apply(paiz, 1, as.list)
      paiz = lapply(paiz, unlist)
      h = l2b_i[,rh100]
      pai = l2b_i[,pai]
      vertprofile = mapply(.pai_vertprofile, paiz, h, pai, SIMPLIFY = F)
      l2b_i$hz = lapply(vertprofile, function(x){unlist(x[,1])})
      l2b_i$paiz = lapply(vertprofile, function(x){unlist(x[,2])})
      gedi[[i]] = l2b_i
    }
  }
  gedi = bind_rows(gedi)
  return(gedi)
}
  
  
  
  
  
  
  
  
  