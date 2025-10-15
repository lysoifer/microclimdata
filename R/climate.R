#' @title Downloads NCEP-NCAR Reanalysis climate data
#'
#' @description This function downloads six-hourly historic climate data from the NCEP-NCAR reanalysis project
#' @param r a SpatRaster giving a template area for which to download data
#' @param tme a POSIXlt object giving the time period for which to download data
#' @return a list of the following:
#' \describe{
#'  \item{tme}{POSIXlt object of period covered by download - typically longer than `tme`}
#'  \item{clim data}{a list of wrapped multi=layer SpatRasters of each climate variable}
#'  \item{tmin.2m}{Minimum temperature (K)}
#'  \item{tmax.2m}{Maximum temperature (K)}
#'  \item{shum.2m}{Specific humidity kg/kg}
#'  \item{pres.sfc}{Sea-level pressure Pa}
#'  \item{dswrf.sfc}{Downward shortwave radiation (W/m^2)}
#'  \item{dlwrf.sfc}{Downward longwave radiation (W/m^2)}
#'  \item{uwnd.10m}{u vector of wind speed at 10m (m/s)}
#'  \item{vwnd.10m}{v vector of wind speed at 10m (m/s)}
#'  \item{prate.sfc}{precipitation rate kg/s}
#' }
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimdata, .registration = TRUE
#' @import terra
#' @import RNCEP
#' @export
#' @rdname ncep_download
#' @seealso [ncep_process()]
ncep_download<-function(r, tme) {
  tz<-attr(tme,"tzone")
  if (tz!= "UTC") stop("time zone in tme must be set to UTC")
  vars<-c("tmin.2m","tmax.2m","shum.2m","pres.sfc","dswrf.sfc","dlwrf.sfc","uwnd.10m","vwnd.10m","prate.sfc")
  climdata<-list()
  for (i in 1:9) {
    a<-.get_ncepone(r, tme, vars[i])
    climdata[[i]]<-wrap(a$roo)
    if (i==1) tme2<-a$tme
  }
  names(climdata)<-vars
  return(list(tme=tme2,climdata=climdata))
}
#' @title Downloads ERA5 climate data
#'
#' @description This function downloads hourly historic climate data from the Copernicus
#' climate data store.
#' @param r a SpatRaster giving a template area for which to download data
#' @param tme a POSIXlt object giving the time period for which to download data
#' @param credentials a data,frame of user credentials (see vignette)
#' @param file_prefix character prefix to attach to the downloaded .nc files when downloaded.
#' The year, month, and file extension is automatically added.
#' @param pathout directory to which data are saved.
#' @param clean if TRUE, removes zip files and only saves nc files
#' @return a list of the request details sent to climateData Store for all microclimate
#' relevant climate variables. Used by `era5_process` to process data. Files a downlaoded
#' to the `pathout` directory.
#' @seealso [era5_process()]
#' @import terra
#' @export
era5_download<-function(r, tme, credentials, file_prefix, pathout, clean = T) {
  # set credentials
  sel <- which(credentials$Site == "CDS")
  uid <- credentials$username[sel]
  cds_access_token <- credentials$password[sel]
  ecmwfr::wf_set_key(user = uid,
                     key = cds_access_token)
  # build request
  e<-ext(r)
  ro<-rast(e)
  crs(ro)<-crs(r)
  rll<-project(ro,"EPSG:4326")
  ell<-ext(rll)
  xmn<-floor(ell$xmin*4)/4
  xmx<-ceiling(ell$xmax*4)/4
  ymn<-floor(ell$ymin*4)/4
  ymx<-ceiling(ell$ymax*4)/4
  req <- mcera5::build_era5_request(xmin = xmn, xmax = xmx,
                                    ymin = ymn, ymax = ymx,
                                    start_time = tme[1],
                                    end_time = tme[length(tme)],
                                    outfile_name = file_prefix)
  # check which downloads already exist
  keep<-rep(TRUE,length(req))
  for (i in 1:length(req)) {
    fi<-paste0(pathout,req[[i]]$target)
    fi<-gsub("zip", "nc", fi) # test if .nc file exists
    if (file.exists(fi)) keep[i]<-FALSE
  }
  s<-which(keep)
  req2<-req[s]
  # download data
  dir.create(pathout,showWarnings=FALSE)
  if(length(req2) != 0) {
    mcera5::request_era5(request = req2, uid = uid, out_path = pathout)
  }
  
  # remove the zip file if clean == TRUE
  if(clean) {
    fpaths = paste0(pathout,sapply(req, "[[", "target"))
    for(i in fpaths) {
      if(file.exists(fpaths[i])) unlink(fpaths[i])
    }
  }
  
  return(req)
}

#' @title Downloads UK Met Office climate data
#'
#' @description This function downloads daily 1km HadUK-Grid historic temperature
#' and precipitation data for the entire UK.
#' @param tme a POSIXlt object giving the time period for which to download data
#' @param credentials a data,frame of user credentials (see vignette)
#' @param pathout directory to which data are saved.
#' @return a data,frame detailing download details and success of all requested files
#' save in `pathout`.
#' @seealso [haduk_process()], [haduk_blend()]
#' @import curl
#' @export
haduk_download <- function(tme, credentials, pathout) {
  dir.create(pathout,showWarnings=FALSE)
  # get user credentials
  sel<-which(credentials$Site == "CEDA")
  cedausr <- credentials$username[sel]
  cedapwd <- credentials$password[sel]
  cedaprot<-"https://"
  cedaserv<-"dap.ceda.ac.uk"
  if(class(tme)[1]!="POSIXlt") stop("tme not POSIXlt class!!")
  dateseq<-seq(tme[1],tme[length(tme)],'month')
  yrs<-lubridate::year(dateseq)
  mnths<-lubridate::month(dateseq)
  varn<-c("rainfall","tasmax","tasmin")
  # Check time period requested
  tmenow<-as.POSIXlt(Sys.time())
  yr<-tmenow$year+1900
  if (any(yrs > yr)) stop("Function downloads observed data. Cannot be used for future years")
  if (any(yrs == yr)) warning("Data may not yet be available for current year")
  if (any(yrs == (yr-1))) warning("Data may not yet be available for previous year")
  if (any(yrs < 1960)) stop("Temperature data not available prior to 1960")
  # Get date text used in file names
  mtxt<-ifelse(mnths<10,paste0("0",mnths),paste0("",mnths))
  daysofmonth<-c(31,28,31,30,31,30,31,31,30,31,30,31)
  mdays<-daysofmonth[mnths]
  mdays<-ifelse((yrs%%4==0 & mnths==2),29,mdays)
  # Create base url
  urlbase<-paste0("ftp.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.2.0.ceda/1km/")
  fullreport_df<-data.frame()
  for (v in varn){
    urlvar<-file.path(urlbase,v,"day","latest")
    fnames<-paste0(v,"_hadukgrid_uk_1km_day_",yrs,mtxt,"01-",yrs,mtxt,mdays,".nc")
    destfiles<-file.path(pathout,fnames)
    dload_urls<-file.path(urlvar,fnames)
    success_df<-curl::multi_download(urls=dload_urls,
                                     destfiles=destfiles,
                                     userpwd = paste0(cedausr,":",cedapwd) )
    if(any(success_df$success==FALSE)) {
      print(paste('UNSUCCESSFUL file downloads:',fnames[which(!success_df$success)]))
    } else print('All downloads successful')
    if(nrow(fullreport_df)==0) {
      fullreport_df<-success_df
    } else fullreport_df<-rbind(fullreport_df,success_df)
  }
  return(fullreport_df)
}
#' @title Processes NCEP-NCAR Reanalysis climate data
#'
#' @description This function processes six-hourly climate data from the NCEp-NCAR
#' reanalysis project for use with `microgrid` or `micropoint`
#' @param ncepdata list of climate variables as returned by `ncep_download` - data
#' downloaded if not supplied
#' @param r a SpatRaster object as a template or downloading data or for resampling
#' and masking outputs if specified.
#' @param tme POSIXlt object for which processed data are required. Must overlap
#' with that supplied in `ncepdata`.
#' @param out one of `grid` or `point` specifying whether to return a list of
#' arrays of climate data or a data.frame for a single point location.
#' @param lat latitude (decimal degrees) - alternative to supplying a template
#' raster if data are required for a point location.
#' @param lon longitude (decimal degrees) - alternative to supplying a template
#' raster if data are required for a point location.
#' @param resampleout optional logical indicating whether to resample and mask outputs
#' to match `r`.
#' @return if `out = "point"` then a a data.frame of hourly weather comprising
#' the following columns:
#' \describe{
#'  \item{obs_time}{POSIXlt object of dates and times}
#'  \item{temp}{Temperature (deg C)}
#'  \item{relhum}{Relative humidity (percentage)}
#'  \item{pres}{Sea-level pressure (kPa)}
#'  \item{swdown}{Downward shortwave radiation on the horizontal (W/m^2)}
#'  \item{difrad}{Downward diffuse radiation (W/m^2)}
#'  \item{lwdown}{Downward longwave radiation W/m^2}
#'  \item{windspeed}{wind speed at 2m (m/s)}
#'  \item{winddir}{wind direction (decimal degrees)}
#'  \item{precip}{precipitation mm}
#' }
#' @return If `out = "grid"`a list of wrapped multilayer SpatRasters of each
#' of the above climate variables (i.e. `tme` not returned - same as input).
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimdata, .registration = TRUE
#' @import terra
#' @export
#' @rdname ncep_process
#' @seealso [ncep_download()]
ncep_process<-function(ncepdata = NA, r, tme, out = "grid", lat = NA, long = NA, resampleout = FALSE) {
  if (class(r) == "logical") {
    r<-rast(matrix(0,ncol=3,nrow=3))
    e<-ext(lon-0.1,lon+0.1,lat-0.1,lon+0.1)
    ext(r)<-e
    crs(r)<-"EPSG:4326"
  }
  if (class(ncepdata) == "logical") ncepdata<-ncep_download(r, tme)
  # These variables are 6 hour hindcasts from the reference time
  # ~~ tmin.2m, tmax.2m
  # These variables are forecasts valid 6 hours after the reference time
  # ~~ "shum.2m","pres.sfc","uwnd.10m","vwnd.10m"
  # These variables are 6 hour averages starting at the reference time
  # ~~ dswrf.sfc", "dlwrf.sfc", prate.sfc
  # ** (1) Get lats and longs and times of array of array
  tmeb<-ncepdata$tme
  climdata<-ncepdata$climdata
  ro<-rast(climdata$tmin.2m)[[1]]
  crs(ro)<-"EPSG:4326"
  if (out == "grid") {
    ll<-.latslonsfromr(ro)
    tc<-rh<-pk<-swr<-dfr<-lwr<-ws<-wdi<-ph<-array(NA,c(dim(ll$lats),length(tme)))
    for (i in 1:dim(ll$lats)[1]) {
      for (j in 1:dim(ll$lons)[2]) {
        ha<-.hourlyall(as.array(rast(climdata$tmin.2m))[i,j,],
                       as.array(rast(climdata$tmax.2m))[i,j,],
                       as.array(rast(climdata$pres.sfc))[i,j,],
                       as.array(rast(climdata$shum.2m))[i,j,],
                       as.array(rast(climdata$dswrf.sfc))[i,j,],
                       as.array(rast(climdata$dlwrf.sfc))[i,j,],
                       as.array(rast(climdata$uwnd.10m))[i,j,],
                       as.array(rast(climdata$vwnd.10m))[i,j,],
                       as.array(rast(climdata$prate.sfc))[i,j,],
                       tmeb,tme,ll$lats[i,j],ll$lons[i,j])
        tc[i,j,]<-ha$temp
        rh[i,j,]<-ha$relhum
        pk[i,j,]<-ha$pres
        swr[i,j,]<-ha$swdown
        dfr[i,j,]<-ha$difrad
        lwr[i,j,]<-ha$lwdown
        ws[i,j,]<-ha$windspeed
        wdi[i,j,]<-ha$winddir
        ph[i,j,]<-ha$prec
      }
    }
    climarrayr<-list()
    climarrayr[[1]]<-.rast(tc,ro)
    climarrayr[[2]]<-.rast(rh,ro)
    climarrayr[[3]]<-.rast(pk,ro)
    climarrayr[[4]]<-.rast(swr,ro)
    climarrayr[[5]]<-.rast(dfr,ro)
    climarrayr[[6]]<-.rast(lwr,ro)
    climarrayr[[7]]<-.rast(ws,ro)
    climarrayr[[8]]<-.rast(wdi,ro)
    climarrayr[[9]]<-.rast(ph,ro)
    if (resampleout) {
      for (i in 1:9) {
        if (crs(ro) != crs(r)) {
          climarrayr[[i]]<-project(climarrayr[[i]],r)
        } else {
          climarrayr[[i]]<-resample(climarrayr[[i]],r)
        }
        climarrayr[[i]]<-mask(climarrayr[[i]],r)
      }
    }
    for (i in 1:9) climarrayr[[i]]<-wrap(climarrayr[[i]])
    names(climarrayr)<-c("temp","relhum","pres","swdown","difrad","lwdown","windspeed","winddir","precip")
  } else {
    ll<-.latlongfromrast(r)
    ha<-.hourlyall(.extractvar(climdata$tmin.2m,ll),
                   .extractvar(climdata$tmax.2m,ll),
                   .extractvar(climdata$pres.sfc,ll),
                   .extractvar(climdata$shum.2m,ll),
                   .extractvar(climdata$dswrf.sfc,ll),
                   .extractvar(climdata$dlwrf.sfc,ll),
                   .extractvar(climdata$uwnd.10m,ll),
                   .extractvar(climdata$vwnd.10m,ll),
                   .extractvar(climdata$prate.sfc,ll),
                   tmeb,tme,ll$lat,ll$long)
    ha<-as.data.frame(ha)
    obs_time<-tme
    climarrayr<-cbind(obs_time,ha)
  }
  return(climarrayr)
}
#' @title Processes ERA5 climate data
#'
#' @description This function processes hourly ERA5 climate data for use with `microgrid`
#' or `micropoint`
#' @param req a list of the request details sent to climateData Store as returned
#' by `era5_download`. If set to NA, all files in directory `pathin` are read in.
#' @param pathin directory in which files downloaded by `era5_download` are saved.
#' @param r a SpatRaster object used as a template for cropping and reprojecting the
#' outputs and resampling and masking outputs if specified.
#' @param out one of `grid` or `point` specifying whether to return a list of
#' arrays of climate data or a data.frame for a single point location.
#' @param lat latitude (decimal degrees) - alternative to supplying a template
#' raster if data are required for a point location.
#' @param lon longitude (decimal degrees) - alternative to supplying a template
#' raster if data are required for a point location.
#' @param resampleout optional logical indicating whether to resample and mask outputs
#' to match `r`.
#' @return if `out = "point"` then a a data.frame of hourly weather comprising
#' the following columns:
#' \describe{
#'  \item{obs_time}{POSIXlt object of dates and times}
#'  \item{temp}{Temperature (deg C)}
#'  \item{relhum}{Relative humidity (percentage)}
#'  \item{pres}{Sea-level pressure (kPa)}
#'  \item{swdown}{Downward shortwave radiation on the horizontal (W/m^2)}
#'  \item{difrad}{Downward diffuse radiation (W/m^2)}
#'  \item{lwdown}{Downward longwave radiation W/m^2}
#'  \item{windspeed}{wind speed at 2m (m/s)}
#'  \item{winddir}{wind direction (decimal degrees)}
#'  \item{precip}{precipitation mm}
#' }
#' @return If `out = "grid"`a list of wrapped SpatRasters of the above climate variables.
#' @import terra
#' @export
#' @rdname era5_process
#' @seealso [era5_download()]
era5_process <- function(req, pathin, r, tme, out = "grid", lat = NA, long = NA, resampleout = FALSE) {
  if (class(req) == "logical") {
    files<-list.files(pathin)
    n<-length(files)
  } else {
    n<-length(req)
    files <- ""
    for (i in 1:n) files[i]<-req[[i]]$target
  }
  # Get raster template for nc file
  files = gsub("zip", "nc", files)
  nc<-paste0(pathin,files[1])
  # Grid process
  climdata<-.extract_clima(nc, r, resampleout)
  if (n > 1) {
    for (i in 2:n) {
      nc<-paste0(pathin,files[i])
      climone<-.extract_clima(nc, r, resampleout)
      for (j in 1:9)  climdata[[j]]<-c(climdata[[j]], climone[[j]])
    }
  }
  tc<-climdata[[1]]
  tmenc<-as.POSIXlt(time(tc),tz="UTC")
  if (length(tmenc) > length(tme)) {
    s<-which(tmenc >= tme[1] & tmenc <= tme[length(tme)])
    for (i in 1:9) {
      ro<-climdata[[i]]
      ro<-ro[[s]]
      climdata[[i]]<-ro
    }
  }
  if (out == "point") {
    if (class(r) != "logical") {
      ll<-.latlongfromrast(r)
      lat <- ll$lat
      long <- ll$long
    }
    xy<-data.frame(x=long,y=lat)
    m<-matrix(0,ncol=10,nrow=length(tme))
    for (i in 1:9) {
      v<-as.numeric(extract(climdata[[i]],xy))
      m[,i+1]<-v[-1]
    }
    climdata<-data.frame(m)
    names(climdata)<-c("obs_time","temp","relhum","pres","swdown","difrad","lwdown","windspeed","winddir","precip")
    climdata$obs_time<-tme
  } else {
    for (i in 1:9)  climdata[[i]]<-wrap(climdata[[i]])
  }
  return(climdata)
}
#' @title Process UK Met Office climate data
#'
#' @description This function processes daily 1km HadUK-Grid historic temperature
#' @param r a SpatRaster object used as a template region for which to extract data
#' from downloaded files.
#' @param tme a POSIXlt object giving the time period for which to download data
#' @param pathin directory to which data are saved using `haduk_download`.
#' @return a list of wrapped SpatRasters of the following:
#' \describe{
#'  \item{precipitation}{daily precipitation (mm)}
#'  \item{tasmax}{Daily maximum temperature (deg C)}
#'  \item{tasmin}{Daily minimum temperature (deg C)}
#' }
#' @seealso [haduk_download()], [haduk_blend()]
#' @import terra
#' @export
haduk_process<-function(r, tme, pathin) {
  # sort out times
  if(class(tme)[1]!="POSIXlt") stop("tme not POSIXlt class!!")
  dateseq<-seq(tme[1],tme[length(tme)],'month')
  yrs<-lubridate::year(dateseq)
  mnths<-lubridate::month(dateseq)
  # Get date text used in file names
  mtxt<-ifelse(mnths<10,paste0("0",mnths),paste0("",mnths))
  daysofmonth<-c(31,28,31,30,31,30,31,31,30,31,30,31)
  mdays<-daysofmonth[mnths]
  mdays<-ifelse((yrs%%4==0 & mnths==2),29,mdays)
  # sort out extents
  ro<-rast(ext(r))
  crs(ro)<-crs(r)
  # create blank OSGB rast
  rte<-rast(matrix(0,ncol=10,nrow=10))
  crs(rte)<-"EPSG:27700"
  if (crs(ro) != crs(rte)) {
    rosgb<-project(ro,"EPSG:27700")
  } else rosgb<-ro
  varn<-c("rainfall","tasmax","tasmin")
  climdata<-list()
  ii<-1
  for (v in varn) {
    fnames<-paste0(v,"_hadukgrid_uk_1km_day_",yrs,mtxt,"01-",yrs,mtxt,mdays,".nc")
    destfiles<-file.path(pathin,fnames)
    rma<-rast(destfiles[1])
    rma<-crop(rma,ext(rosgb))
    if (length(destfiles) > 1) {
      for (i in 2:length(destfiles)) {
        ri<-rast(destfiles[i])
        ri<-crop(ri,ext(rosgb))
        rma<-c(rma,ri)
      }
    }
    climdata[[ii]]<-wrap(rma)
    ii<-ii+1
  }
  names(climdata)<-c("precipitation","tasmax","tasmin")
  return(climdata)
}
#' @title Combine coarse-resolution and UK Met Office data
#' @description This function blends coarse resolution and daily 1km HadUK-Grid
#' historic temperature  and precipitation data to provide hourly 1km climate data
#' @param coarsedata coarse-resolution hourly climate data as returned by `ncep_process`
#' or `era5_process`.
#' @param hadukdata 1km daily precipitation and temperature data as returned by
#' `haduk_process`
#' @return a list comprising `climarrayr` - wrapped SpatRasters of 1km hourly
#' climate data and `tme` - a POSIXlt object of dates and times corresponding to
#' data in `climarrayr`. The dataset `climarrayr` is itself a list of wrapped
#' SpatRasters of the following variables:
#' \describe{
#'  \item{temp}{Temperature (deg C)}
#'  \item{relhum}{Relative humidity (percentage)}
#'  \item{pres}{Sea-level pressure (kPa)}
#'  \item{swdown}{Downward shortwave radiation on the horizontal (W/m^2)}
#'  \item{difrad}{Downward diffuse radiation (W/m^2)}
#'  \item{lwdown}{Downward longwave radiation W/m^2}
#'  \item{windspeed}{wind speed at 2m (m/s)}
#'  \item{winddir}{wind direction (decimal degrees)}
#'  \item{precip}{precipitation mm}
#' }
#' @seealso [haduk_download()], [haduk_process()]
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimdata, .registration = TRUE
#' @export
haduk_blend<-function(coarsedata, hadukdata) {
  # check temporal
  nh<-dim(rast(coarsedata$temp))[3]
  nd<-dim(rast(hadukdata$tasmin))[3]
  if (nh != 24*nd) stop("coarsedata does not match hadukdata temporally")
  # check spatial
  rre<-rast(coarsedata$temp)[[1]]
  rte<-rast(hadukdata$tasmin)[[1]]
  if (crs(rre) != crs(rte)) {
    r1<-project(rre,crs(rte))
  } else r1<-rre
  ist<-terra::intersect(ext(r1),ext(rte))
  if (class(ist) == "logical")  stop("spatial extent of coarsedata does not overlap with hadukdata")
  # Temperature
  tasmin<-rast(hadukdata$tasmin)
  tasmax<-rast(hadukdata$tasmax)
  temp<-rast(coarsedata$temp)
  # reproject / resample
  if (crs(temp) != crs(tasmin)) {
    temp<-project(temp,rte)
  } else {
    temp<-resample(temp,rte)
  }
  temp<-mask(temp,tasmin[[1]])
  e<-ext(temp)
  tasmin<-crop(tasmin,e)
  tasmax<-crop(tasmax,e)
  tasmin<-mask(tasmin,temp[[1]])
  tasmax<-mask(tasmax,temp[[1]])
  temp<-blendtempCpp(.is(tasmin),.is(tasmax),.is(temp))
  temp<-.rast(temp,rte)
  # Precipitation
  precd<-rast(hadukdata$precipitation)
  prech<-rast(coarsedata$prec)
  # reproject / resample
  if (crs(prech) != crs(precd)) {
    prech<-project(prech,precd[[1]])
  } else {
    prech<-resample(prech,precd[[1]])
  }
  prech<-mask(prech,precd[[1]])
  e<-ext(prech)
  precd<-crop(precd,e)
  precd<-mask(precd,prech[[1]])
  prech<-blendprecCpp(.is(precd),.is(prech))
  prech<-.rast(prech,precd[[1]])
  # Relative humidity
  tc<-.is(rast(coarsedata$temp))
  rh<-.is(rast(coarsedata$relhum))
  ea<-.rast(satvapCpp(tc)*rh/100,rre)
  if (crs(ea) != crs(rte)) {
    ea<-project(ea,rte)
  } else {
    ea<-resample(ea,rte)
  }
  rh<-(.is(ea)/satvapCpp(.is(temp)))*100
  rh[rh>100]<-100
  rh[rh<0]<-0
  rh<-.rast(rh,rte)
  # U and V wind
  u2<-.is(coarsedata$windspeed)*cos(.is(coarsedata$winddir)*pi/180)
  v2<-.is(coarsedata$windspeed)*sin(.is(coarsedata$winddir)*pi/180)
  u2<-.rast(u2,rre)
  v2<-.rast(v2,rre)
  if (crs(u2) != crs(rte)) {
    u2<-project(u2,rte)
    v2<-project(v2,rte)
  } else {
    u2<-resample(u2,rte)
    v2<-resample(v2,rte)
  }
  ws<-sqrt(.is(u2)^2+.is(v2)^2)
  wd<-(atan2(.is(u2),.is(v2))*180/pi+180)%%360
  ws<-.rast(ws,rte)
  wd<-.rast(wd,rte)
  ws<-mask(ws,rte)
  wd<-mask(wd,rte)
  # Shortwave radiation
  tme<-as.POSIXlt(time(rast(coarsedata$temp)),tz="UTC")
  year<-tme$year+1900
  month<-tme$mon+1
  day<-tme$mday
  hour<-tme$hour
  ll<-.latslonsfromr(rre)
  csr<-clearskya(year,month,day,hour,ll$lats, ll$lons)
  csf<-rast(coarsedata$swdown)/.rast(csr,rre)
  if (crs(csf) != crs(rte)) {
    csf<-project(csf,rte)
  } else {
    csf<-resample(csf,rte)
  }
  csf[is.na(csf)]<-0
  csf[csf>1]<-1
  ll<-.latslonsfromr(rte)
  csr<-clearskya(year,month,day,hour,ll$lats,ll$lons)
  swrad<-csf*.rast(csr,rte)
  swrad<-mask(swrad,rte)
  # Diffuse rasdiation
  dsf<-rast(coarsedata$difrad)/rast(coarsedata$swdown)
  dsf[is.na(dsf)]<-0.5
  if (crs(dsf) != crs(rte)) {
    dsf<-project(dsf,rte)
  } else {
    dsf<-resample(dsf,rte)
  }
  difrad<-dsf*swrad
  difrad<-mask(difrad,rte)
  s<-which(.is(difrad) > .is(swrad))
  difrad[s]<-swrad[s]
  # lwrad
  skyem<-rast(coarsedata$lwdown)/(5.67*10^-8*(rast(coarsedata$temp)+273.15)^4)
  if (crs(skyem) != crs(rte)) {
    skyem<-project(skyem,rte)
  } else {
    skyem<-resample(skyem,rte)
  }
  skyem[skyem>1]<-1
  skyem[skyem<0]<-0
  lwrad<-skyem*5.67*10^-8*(temp+273.15)^4
  # pressure
  pk<-rast(coarsedata$pres)
  if (crs(pk) != crs(rte)) {
    pk<-project(pk,rte)
  } else {
    pk<-resample(pk,rte)
  }
  pk<-mask(pk,rte)
  # Create output
  climarrayr<-list()
  climarrayr[[1]]<-temp
  climarrayr[[2]]<-rh
  climarrayr[[3]]<-pk
  climarrayr[[4]]<-swrad
  climarrayr[[5]]<-difrad
  climarrayr[[6]]<-lwrad
  climarrayr[[7]]<-ws
  climarrayr[[8]]<-wd
  climarrayr[[9]]<-prech
  # add times and wrap
  for (i in 1:9) {
    time(climarrayr[[i]])<-as.POSIXct(tme)
    climarrayr[[i]]<-wrap(climarrayr[[i]])
  }
  names(climarrayr)<-c("temp","relhum","pres","swdown","difrad","lwdown","windspeed","winddir","precip")
  return(climarrayr)
}
#' @title Fills missing coastal grid cells in UK Met Office data
#' @description This function infills missing coastal grid cells of the 1 km HadUK-Grid
#' @param climdata a list of climate data as returned by `haduk_blend`
#' @param sst a multi-layer SpatRast of monthly sst temperature data corresponding
#' to `climdata`
#' @param wgt an optional weighting factor between 0 and 1 indicating the sensitivity
#' of coastal temperatures to sea-surface temperatures.
#' @return a list comprising `climarrayr` - wrapped SpatRasters of 1km hourly
#' climate data and `tme` - a POSIXlt object of dates and times corresponding to
#' data in `climarrayr`. The dataset `climarrayr` is itself a list of wrapped
#' SpatRasters of the following variables:
#' \describe{
#'  \item{temp}{Temperature (deg C)}
#'  \item{relhum}{Relative humidity (percentage)}
#'  \item{pres}{Sea-level pressure (kPa)}
#'  \item{swdown}{Downward shortwave radiation on the horizontal (W/m^2)}
#'  \item{difrad}{Downward diffuse radiation (W/m^2)}
#'  \item{lwdown}{Downward longwave radiation W/m^2}
#'  \item{windspeed}{wind speed at 2m (m/s)}
#'  \item{winddir}{wind direction (decimal degrees)}
#'  \item{precip}{precipitation mm}
#' }
#' @seealso [haduk_blend()]
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimdata, .registration = TRUE
#' @export
haduk_coastfill<-function(climdata, sst, wgt = 0.5) {
  # ==================================================================== #
  # ~~~~~~~~~~~~~~~~ get lsrin data  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  # ==================================================================== #
  lsrr<-rast(ukcoastexposure[[1]])
  for (i in 2:8) lsrr<-c(lsrr,rast(ukcoastexposure[[i]]))
  # ==================================================================== #
  # ~~~~~~~~~~~ Clean sst data and create coastal layer  ~~~~~~~~~~~~~~~ #
  # ==================================================================== #
  tc<-rast(climdata$temp)
  rte<-tc[[1]]
  sst<-.expand_sst(sst)
  if (crs(r) != crs(sst)) {
    sst<-project(sst,r)
  } else sst<-resample(sst,r)
  # crop lsr
  lsrr<-.is(crop(lsrr,ext(rte)))
  lsm<-applymeanCpp(lsrr)
  lsea<-lsm*0+1
  # ==================================================================== #
  # ~~~~~~~ Sort out non-temperature dependent variables ~~~~~~~~~~~~~~~ #
  # ==================================================================== #
  # pressure
  pk<-.is(rast(climdata$pres))
  pk<-.rast(na_fill_array(pk,.is(lsea),2),rte)
  # shortwave radiation
  tme<-as.POSIXlt(time(rast(climdata$temp)),tz="UTC")
  year<-tme$year+1900
  month<-tme$mon+1
  day<-tme$mday
  hr<-tme$hour+tme$min/60+tme$sec/3600
  ll<-.latslonsfromr(rte)
  csr<-clearskya(year,month,day,hr,ll$lats,ll$lons)
  csf<-.is(climdata$swdown)/csr
  csf[csr==0]<-0
  csf<-mask(.rast(csf,rte),rte)
  csf<-na_fill_array(.is(csf),.is(lsea),2)
  swdown<-.rast(csf*csr,rte)
  # diffuse radiation
  difr<-.is(climdata$difrad)/.is(climdata$swdown)
  difr[csr==0]<-0
  difr<-mask(.rast(difr,rte),rte)
  difr<-na_fill_array(.is(difr),.is(lsea),2)
  difr<-.rast(difr*.is(swdown),rte)
  # precipitation
  prec<-.is(rast(climdata$precip))
  prec<-.rast(na_fill_array(prec,.is(lsea),2),rte)
  prec[prec<0]<-0
  # ==================================================================== #
  # ~~~~~~~~~~~~~~~~~~~~~~ Calculate modal wind ~~~~~~~~~~~~~~~~~~~~~~~~ #
  # ==================================================================== #
  # Calculate modal wind
  wd<-rast(climdata$winddir)
  mx<-min(dim(wd)[1:2])
  if (mx > 5) wd<-aggregate(wd,5,mean,na.rm=TRUE)
  wdir<-modalwinddir(.is(wd))
  # Calculate array land-sea ratios for every hour
  i<-round(wdir/45)%%8
  lsr<-lsrr[,,i+1]
  # fill wind speed and direction
  u2<-.is(climdata$windspeed)*cos(.is(climdata$winddir)*pi/180)
  v2<-.is(climdata$windspeed)*sin(.is(climdata$winddir)*pi/180)
  u2<-na_fill_array(u2,.is(lsea),2)
  v2<-na_fill_array(v2,.is(lsea),2)
  ws<-.rast(sqrt(u2^2+v2^2),rte)
  wd<-.rast((atan2(u2,v2)*180/pi+180)%%360,rte)
  # ==================================================================== #
  # ~~~~~~~~~~~~~~~ Calculate coastal temperature effect ~~~~~~~~~~~~~~~ #
  # ==================================================================== #
  # expand tc
  tc<-.rast(na_fill_array(.is(tc),.is(lsea),2),rte)
  # derive power scaling coefficient from wind speed
  p2<-0.10420*sqrt(.is(ws))-0.47852
  # calculate logit lsr and lsm
  llsr<-log(lsr/(1-lsr))
  llsm<-log(lsm/(1-lsm))
  llsm<-mattoarray(llsm,dim(llsr)[3])
  # Calculate mins and maxes
  s<-which(is.na(lsr) == FALSE & lsr > 0 & lsr < 1)
  llsr[lsr==0]<-log(min(lsr[s])/(1-min(lsr[s])))
  llsr[lsr==1]<-log(max(lsr[s])/(1-max(lsr[s])))
  s<-which(is.na(lsm) == FALSE & lsm > 0 & lsm < 1)
  llsm[lsm==0]<-log(min(lsm[s])/(1-min(lsm[s])))
  llsm[lsm==1]<-log(max(lsm[s])/(1-max(lsm[s])))
  # predict logit swgt
  lswgt<- -0.1095761+p2*(llsr+3.401197)-0.1553487*llsm
  swgt<-.rast(1/(1+exp(-lswgt)),rte)
  swgt<-wgt*swgt
  tcp<-swgt*sst+(1-swgt)*tc
  # calculate aggregation factor
  a1<-.is(rast(climdata$temp))
  a2<-.is(tcp)
  s<-which(is.na(a1) & is.na(a2) == FALSE)
  a1[s]<-a2[s]
  tc<-.rast(a1,rte)
  # ==================================================================== #
  # ~~~~~~~~~~~~~~~ Calculate temperature-dependent variables ~~~~~~~~~~ #
  # ==================================================================== #
  # relative humidity
  ea<-satvapCpp(.is(rast(climdata$temp)))*.is(rast(climdata$relhum))/100
  ea<-na_fill_array(ea,.is(lsea),2)
  rh<-(ea/satvapCpp(.is(tc)))*100
  rh[rh>100]<-100
  # longwave radiation
  lwup<-5.67*10^-8*(.is(rast(climdata$temp))+273.15)^4
  skyem<-.is(rast(climdata$lwdown))/lwup
  skyem<-na_fill_array(.is(skyem),.is(lsea),2)
  lwdown<-.rast(skyem*5.67*10^-8*(.is(tc)+273.15)^4,rte)
  # turn into list
  climarrayr<-list()
  climarrayr$temp<-tc
  climarrayr$relhum<-.rast(rh,rte)
  climarrayr$pres<-pk
  climarrayr$swdown<-swdown
  climarrayr$difrad<-difr
  climarrayr$lwdown<-lwdown
  climarrayr$windspeed<-ws
  climarrayr$winddir<-wd
  climarrayr$precip<-prec
  # add times and wrap
  for (i in 1:9) {
    time(climarrayr[[i]])<-as.POSIXct(tme)
    climarrayr[[i]]<-wrap(climarrayr[[i]])
  }
  return(climarrayr)
}
#' @title Downloads daily sea-surface temperature data
#' @description This function downloads daily 0.25 degree grid resolution sea-surface
#' temperature data from the NOAA National Centers for Environmental Information
#' @param r a SpatRaster giving a template area for which to download data
#' @param tme a POSIXlt object giving the time period for which to download data
#' @return a multi-layer SpatRaster of daily sea-surface temperature data
#' @seealso [haduk_blend()]
#' @import terra
#' @import rerddap
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimdata, .registration = TRUE
#' @export
sst_download<-function(r, tme, resampleout = "FALSE", nafill = "FALSE") {
  ## Define lat and long range
  rll<-rast(r)
  ext(rll)<-ext(r)
  crs(rll)<-crs(r)
  rll<-project(rll,"EPSG:4326")
  e<-ext(rll)
  latitude_range<-c(e$ymin,e$ymax)
  longitude_range<-c(e$xmin,e$xmax)
  # define time range
  st_time <- format(tme[1], format = "%Y-%m-%d")
  ed_time <- format(tme[length(tme)], format = "%Y-%m-%d")
  sst_data <- griddap(
    datasetx = "ncdcOisst21Agg_LonPM180",
    latitude = latitude_range,
    longitude = longitude_range,
    time = c(st_time, ed_time),
    fields = "sst",
    fmt = "nc"
  )
  sst<-sst_data$data
  tms<- unique(sst$time)
  sst$zlev<-NULL
  s<-which(sst$time==tms[1])
  ssto<-sst[s,]
  ssto<-data.frame(x=ssto$longitude,y=ssto$latitude,z=ssto$sst)
  sstm<-rast(ssto, type="xyz")
  for (i in 2:length(tms)) {
    s<-which(sst$time==tms[i])
    ssto<-sst[s,]
    ssto<-data.frame(x=ssto$longitude,y=ssto$latitude,z=ssto$sst)
    ssto<-rast(ssto, type="xyz")
    sstm<-c(sstm,ssto)
  }
  names(sstm)<-tms
  crs(sstm)<-"EPSG:4326"
  if (nafill) {
    rte<-sstm[[1]]*0
    rte[is.na(rte)]<-0
    sstm<-.rast(na_fill_array(.is(sstm),.is(rte)),rte)
  }
  if (resampleout) {
    if(crs(r) != crs(sstm)) {
      sstm<-project(sstm,r)
    } else sstm<-resample(sstm,r)
  }
  return(sstm)
}
#' @title extracts climate data for a point location
#' @description This function extracts climate data from a climate arrays for any
#' given point location returning a data,frame as used by the `micropoint` package
#' @param climdata a list of SpatRasters of climate data as returned by haduk_coastfill or
#' any other of the climate data process functions.
#' @param x easting or longitude
#' @param y northing or latitude
#' @param aslatlong optional logical indicating whether `x` and `y` are longitude
#' and latitude (TRUE) or the easting and northing in the same coordinate reference
#' system as `climdata`
#' @import terra
#' @export
climpoint_extract<-function(climdata, x, y, aslatlong = FALSE) {
  r<-rast(climdata[[1]])[[1]]
  xy <- data.frame(x = x, y = y)
  if (aslatlong) {
    xy <- sf::st_as_sf(xy, coords = c("x", "y"),
                       crs = 4326)
    xy <- sf::st_transform(xy, crs(r))
  } else {
    xy <- sf::st_as_sf(xy, coords = c("x", "y"),
                       crs = crs(r))
  }
  rv<-rast(climdata[[1]])
  dfout<-data.frame(V1=as.POSIXlt(time(rv)))
  for (i in 1:9) {
    rv<-rast(climdata[[i]])
    dfout[,i+1]<-as.numeric(extract(rv, xy))[-1]
  }
  names(dfout)<-c("obs_time",names(climdata))
  return(dfout)
}
#' @title Download UKCP18 climate data, including future climate
#' @description The function downloads specified data from ftp.ceda.ac.uk to a specified directory
#' @param tme POSIXlt object indicating for which time period data are downloaded
#' @param cedatoken ceda access token as provided by \href{https://services-beta.ceda.ac.uk/account/token/?_ga=2.221117848.1874942177.1727205775-405001952.1710779653}{Ceda access token generator}
#' @param pathout directory to which to save data
#' @param collection text string defining UKCP18 collection, either `land-gcm`
#' (60 km grid resolution) or `land-rcm` (12 km grid resolution)
#' @param domain text string defining UKCP18 domain, either 'uk' (United Kingdom)
#' or 'eur' (Europe) or for `land-gcm`, `global` (global extent)
#' @param rcp text string of RCP scenario to be downloaded (either `rcp85` or `rcp26`)
#' @param member vector of strings defining the climate model member to be downloaded. Available members vary between UKCP18 collections.
#' @param silent optional logical indicating whether to print which files are being downloaded.
#' @import curl
#' @export
ukcp18_download<-function(tme, cedatoken, pathout,
                          collection=c('land-gcm','land-rcm'),
                          domain=c('uk','eur','global'),
                          rcp=c('rcp85', 'rcp26'),
                          member=c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15',
                                   '16','17','18','19','20','21','22','23','24','25','26','27','28'),
                          silent = FALSE)
{
  # create list of variables
  vars=c('huss','pr','psl','rls','rss','tasmax','tasmin','uas','vas')
  # create output directory
  dir.create(pathout,showWarnings=FALSE)
  # Check parameters
  collection<-match.arg(collection)
  domain<-match.arg(domain)
  rcp<-match.arg(rcp)
  member<-match.arg(member,several.ok=TRUE)
  # original warnings
  if(collection=='land-rcm' & any(member %in% c('02','03','16','17','18','19','20','21','22','23','24','25','26','27','28'))){
    warning('Invalid member for land-rcm - retricting to only those valid!!')
    member<-member[which(member %in% c('01','04','05','06','07','08','09','10','11','12','13','15'))]
  }
  if(collection=='land-gcm' & !domain %in% c('uk','global')){
    warning('Invalid area for global data - downloading global data!!')
    domain<-'global'
  }
  if(collection=='land-rcm' & !domain %in% c('uk','eur')){
    warning('Invalid area for regional data - downloading Europe data!!')
    domain<-'eur'
  }
  # new warnings
  if (domain=='uk') {
    if(rcp=='rcp85' & any(member %in% c('17','18','22'))){
      warning('Invalid member for land-gcm and UK - retricting to only those valid!!')
      member<-member[which(member %in% c('01','02','03','04','05','06','07','08','09','10',
                                         '11','12','13','14','15','16','19','20','21','23','24','25','26','27','28'))]
    }
    if(rcp=='rcp26' & any(member %in% c('16','17','18','20','21','22','23','24','26','27','28'))){
      warning('Invalid member for land-gcm and UK - retricting to only those valid!!')
      member<-member[which(member %in% c('01','02','03','04','05','06','07','08','09','10','11','12',
                                         '13','14','15','19','25'))]
    }
  }
  if (collection=='land-rcm' & rcp=='rcp26') stop("Invalid rcp for land-gcm")
  # Identify which decades are required
  startdate<-tme[1]
  enddate<-tme[length(tme)]
  decades<-.find_ukcp_decade(collection,startdate,enddate)
  # Get tiled data according to resolution
  if(collection=='land-rcm') reso<-"12km"
  if(collection=='land-gcm') reso<-"60km"
  urlin<-file.path("https://dap.ceda.ac.uk","badc","ukcp18","data",collection,domain,reso,rcp)
  # Download files - loop over model runs and vars - use lapply to download all decades
  output<-c()
  for(run in member){
    for(v in vars){
      fnames<-paste0(v,'_',rcp,'_',collection,'_',domain,'_',reso,'_',run,'_day_',decades,'.nc')
      dload_urls<-file.path(urlin,run,v,'day','latest',fnames)
      destfiles<-file.path(pathout,fnames)
      output<-c(output,destfiles)
      for(n in 1:length(dload_urls)) {
        fi<-destfiles[n]
        if (file.exists(fi) == FALSE) {
          if (silent == FALSE) print(paste('Downloading',run,v,'file',n))
          h <- new_handle(verbose = FALSE)
          handle_setheaders(h, .list=list("Authorization" = paste('Bearer',cedatoken)))
          curl_download(dload_urls[n],destfile=destfiles[n], handle = h)
        }
      }
    }
  }
  # Download orography if requested
  # Checks expected downloaded files exist
  if(any(file.exists(output))==FALSE){
    warning('Download files not found!!!')
    print(output[which(file.exists(output)==FALSE)])
  }
  # Checks downloaded files >100MB
  info<-file.info(output)
  if(any(info$size<100000000)){
    warning('Downloaded files not expected size!!!')
    print(output[which(info$size<100000000)])
  }
  return(output)
}
#' @title Process UKCP18 data converting to daily data
#' @description The function reads in downloaded UKCP18 data, crops out the right area
#' and time slice and conbverts data form 360 days per year to true daily data.
#' @param r SpatRaster giving extent of returned data. If resampleout is set to TRUE,
#' returned data are projected / resampled to match r,
#' @param tme POSIXlt object indicating for which time period data are downloaded
#' @param pathin directory to which data are saved.
#' @param collection text string defining UKCP18 collection, either `land-gcm`
#' (60 km grid resolution) or `land-rcm` (12 km grid resolution)
#' @param domain text string defining UKCP18 domain, either 'uk' (United Kingdom)
#' or 'eur' (Europe) or for `land-gcm`, `global` (global extent)
#' @param rcp text string of RCP scenario to be downloaded (either `rcp85` or `rcp26`)
#' @param member strings defining the climate model member downloaded.
#' @param resampleout optional logical indicating whether to reproject / resample output to match `r`
#' @import terra
#' @export
ukcp18_process<-function(r, tme, pathin,
                         collection=c('land-gcm','land-rcm'),
                         domain=c('uk','eur','global'),
                         rcp=c('rcp85', 'rcp26'),
                         member=c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15',
                                  '16','17','18','19','20','21','22','23','24','25','26','27','28'),
                         resampleout = FALSE) {
  # get list of vars
  vars=c('huss','pr','psl','rls','rss','tasmax','tasmin','uas','vas')
  # check whether anything more than one
  if (length(collection) > 1) stop("collection must be of length one")
  if (length(domain) > 1) stop("domain must be of length one")
  if (length(rcp) > 1) stop("rcp must be of length one")
  if (length(member) > 1) stop("function only works on one member at a time. Please specify single member\n")
  # Identify which decades are required
  startdate<-tme[1]
  enddate<-tme[length(tme)]
  decades<-.find_ukcp_decade(collection,startdate,enddate)
  if(collection=='land-rcm') reso<-"12km"
  if(collection=='land-gcm') reso<-"60km"
  # check whether files exist
  for(v in vars){
    fnames<-paste0(v,'_',rcp,'_',collection,'_',domain,'_',reso,'_',member,'_day_',decades,'.nc')
    destfiles<-file.path(pathin,fnames)
    for (i in 1:length(destfiles)) if (file.exists(destfiles[i]) == FALSE) stop(paste(destfiles[i],"doesn't exist"))
  }
  # Convert projection of r
  rosgb<-project(r,"EPSG:27700")
  e<-ext(rosgb)
  rma<-list()
  j<-1
  if (resampleout) {
    pb<-utils::txtProgressBar(min = 0, max = 18, style = 3)
  } else pb<-utils::txtProgressBar(min = 0, max = 9, style = 3)
  for (v in vars) {
    fnames<-paste0(v,'_',rcp,'_',collection,'_',domain,'_',reso,'_',member,'_day_',decades,'.nc')
    nc<-file.path(pathin,fnames)
    ri<-rast(nc[1],subds=v)
    ri<-crop(ri,e,snap='out')
    # fix times
    ukdates<-.get_ukcp18_dates(nc[1])
    real_dates<-.correct_ukcp_dates(ukdates)
    ri<-.fill_calendar_data(ri,real_dates,testplot=FALSE)
    if (length(nc) > 1) {
      for (i in 2:length(nc)) {
        ri1<-rast(nc[i],subds=v)
        ri1<-crop(ri1,e,snap='out')
        # fix times
        ukdates<-.get_ukcp18_dates(nc[i])
        real_dates<-.correct_ukcp_dates(ukdates)
        ri1<-.fill_calendar_data(ri1,real_dates,testplot=FALSE)
        ri<-c(r1,ri1)
      }
    }
    # subset times
    tme2<-as.POSIXlt(terra::time(ri),tz="UTC")
    s<-which(tme2 >= tme[1] & tme2 <= tme[length(tme)])
    rma[[j]]<-ri[[s]]
    j<-j+1
    utils::setTxtProgressBar(pb,j)
  }
  if (resampleout) {
    for (i in 1:9) {
      if (crs(ro) != crs(r)) {
        rma[[i]]<-project(rma[[i]],r)
      } else {
       rma[[i]]<-resample(rma[[i]],r)
      }
      climarrayr[[i]]<-mask(rma[[i]],r)
      utils::setTxtProgressBar(pb,i+9)
    }
  }
  for (i in 1:9) rma[[i]]<-wrap(rma[[i]])
  names(rma)<-vars
  return(rma)
}



ukcp18_biascorrect<-function(haduk, ukcp_hist, ukcp_futu, albedo = NA, silent = FALSE) {
  r<-rast(haduk$temp)[[1]]
  # ================ Observed hourly to daily =============================== #
  if (silent == FALSE) {
    cat("Converting observed hourly data to daily\n")
    pb<-utils::txtProgressBar(min = 0, max = 8, style = 3)
  }
  # Temperature
  if (silent == FALSE) utils::setTxtProgressBar(pb,0)
  tcc<-temptoday(.is(rast(haduk$temp)))
  tmean<-.rast(tcc$tmean,r)
  if (silent == FALSE) utils::setTxtProgressBar(pb,1)
  dtr<-.rast(tcc$dtr,r)
  if (silent == FALSE) utils::setTxtProgressBar(pb,2)
  # vapour pressure
  ea<-.rast(rhtoday(.is(rast(haduk$relhum)),.is(rast(haduk$temp))),r)
  if (silent == FALSE) utils::setTxtProgressBar(pb,3)
  # atmospheric pressure
  pk<-.rast(applytodaily(.is(rast(haduk$pres)),"mean"),r)
  if (silent == FALSE) utils::setTxtProgressBar(pb,4)
  # swdown
  swdown<-.rast(applytodaily(.is(rast(haduk$swdown)),"mean"),r)
  if (silent == FALSE) utils::setTxtProgressBar(pb,5)
  # sky emissivity
  skyem<-.rast(lwtoday(.is(rast(haduk$lwdown)),.is(rast(haduk$temp))),r)
  if (silent == FALSE) utils::setTxtProgressBar(pb,6)
  # wind speed
  ws<-.rast(applytodaily(.is(rast(haduk$windspeed)),"mean"),r)
  if (silent == FALSE) utils::setTxtProgressBar(pb,7)
  # precipitation
  prec<-.rast(applytodaily(.is(rast(haduk$precip)),"sum"),r)
  if (silent == FALSE) utils::setTxtProgressBar(pb,8)
  # ================ UKCP18 variable conversion =============================== #
  if (silent == FALSE) cat("\nConverting UKCP climate variables\n")
  # temperature
  tcc<-tmnmxtodtr(.is(rast(ukcp_hist$tasmin)), .is(rast(ukcp_hist$tasmax)))
  ukcphist_tmean<-.rast(tcc$tmean,r)
  ukcphist_dtr<-.rast(tcc$dtr,r)
  tcc<-tmnmxtodtr(.is(rast(ukcp_futu$tasmin)), .is(rast(ukcp_futu$tasmax)))
  ukcpfutu_tmean<-.rast(tcc$tmean,r)
  ukcpfutu_dtr<-.rast(tcc$dtr,r)
  # vapour pressure
  ukcphist_ea<-.rast(sphtoea(.is(rast(ukcp_hist$huss)),.is(rast(ukcp_hist$psl))/10),r)
  ukcpfutu_ea<-.rast(sphtoea(.is(rast(ukcp_futu$huss)),.is(rast(ukcp_futu$psl))/10),r)
  # shortwave radiation
  if (class(albedo) == "logical") albedo<-rast(ukcp_hist$rss)*0+0.23
  ukcphist_rsw<-.rast(netshorttodownshort(.is(rast(ukcp_hist$rss)),.is(albedo)),r)
  ukcpfutu_rsw<-.rast(netshorttodownshort(.is(rast(ukcp_futu$rss)),.is(albedo)),r)
  # sky emissivity
  ukcphist_skyem<-.rast(netlwtoskyem(.is(rast(ukcp_hist$rls)),.is(rast(ukcp_hist$tasmax)),.is(rast(ukcp_hist$tasmin))),r)
  ukcpfutu_skyem<-.rast(netlwtoskyem(.is(rast(ukcp_futu$rls)),.is(rast(ukcp_futu$tasmax)),.is(rast(ukcp_futu$tasmin))),r)
  # wind speed
  ukcphist_ws<-.rast(uvtows(.is(rast(ukcp_hist$uas)),.is(rast(ukcp_hist$vas))),r)
  ukcpfutu_ws<-.rast(uvtows(.is(rast(ukcp_futu$uas)),.is(rast(ukcp_futu$vas))),r)
  # ================ UKCP18 bias correction =============================== #
  if (silent == FALSE) cat("\nBias correcting mean daily temperature\n")
  ukcpfutu_tmean<-.biascorrect(tmean, ukcphist_tmean, ukcpfutu_tmean, rangelims, silent)
  if (silent == FALSE) cat("\nBias correcting diurnal temperature range\n")
  ukcpfutu_dtr<-.biascorrect(dtr, ukcphist_dtr, ukcpfutu_dtr, rangelims, silent)
  if (silent == FALSE) cat("\nBias correcting vapour pressure\n")
  ukcpfutu_ea<-.biascorrect(ea, ukcphist_ea, ukcpfutu_ea, rangelims, silent)
  if (silent == FALSE) cat("\nBias correcting atmospheric pressure\n")
  ukcpfutu_pk<-.biascorrect(pk,rast(ukcp_hist$psl)/10,rast(ukcp_futu$psl)/10, rangelims, silent)
  if (silent == FALSE) cat("\nBias correcting shortwave radiation\n")
  ukcpfutu_rsw<-.biascorrect(swdown,ukcphist_rsw,ukcpfutu_rsw, rangelims, silent)
  if (silent == FALSE) cat("\nBias correcting sky emissivity\n")
  ukcpfutu_skyem<-.biascorrect(skyem,ukcphist_skyem,ukcpfutu_skyem, rangelims, silent)
  ukcpfutu_skyem[ukcpfutu_skyem>1]<-1
  ukcpfutu_skyem[ukcpfutu_skyem<0]<-0
  if (silent == FALSE) cat("\nBias correcting wind speed\n")
  ukcpfutu_ws<-.biascorrect(ws,ukcphist_ws,ukcpfutu_ws, rangelims, silent)
  if (silent == FALSE) cat("\nBias correcting precipitation\n")
  ukcpfutu_prec<-.precipcorrect(prec,rast(ukcp_hist$pr),rast(ukcp_futu$pr), rangelims)
  # ================ Expand to hourly ======================================= #
  if (silent == FALSE) {
    cat("\nExpanding bias-corrected daily values to hourly\n")
    pb<-utils::txtProgressBar(min = 0, max = 15, style = 3)
  }
  # Day stuff
  tmed<-as.POSIXlt(time(rast(ukcp_futu$tasmax)))
  year<-tmed$year+1900
  month<-tmed$mon+1
  day<-tmed$mday
  ll<-.latslonsfromr(r)
  # Temperature
  if (silent == FALSE) utils::setTxtProgressBar(pb,1)
  temp<-.rast(tempha(.is(ukcpfutu_tmean), .is(ukcpfutu_dtr),
                     year, month, day, ll$lats, ll$lons),r)
  # Pressure
  if (silent == FALSE) utils::setTxtProgressBar(pb,2)
  pk<-.rast(splina(.is(ukcpfutu_pk)),r)
  # Relative humidity
  if (silent == FALSE) utils::setTxtProgressBar(pb,3)
  rh<-.rast(relhuma(.is(ukcpfutu_ea),.is(temp)),r)
  # SW down
  # ~ clear-sky rad  as function of latitude
  mnlat<-min(ll$lats)-1
  mxlat<-max(ll$lats)+1
  mnlon<-min(ll$lons)-1
  mxlon<-max(ll$lons)+1
  lats<-seq(mxlat,mnlat,length.out=dim(r)[[1]]+2)
  csr<-dailyclm(year, month, day, lats)
  csr<-rast(dailycla(csr, dim(r)[[2]]+2))
  ext(csr)<-ext(mnlon,mxlon,mnlat,mxlat)
  crs(csr)<-"EPSG:4326"
  csr<-project(csr,r)
  if (silent == FALSE) utils::setTxtProgressBar(pb,4)
  # ~ calculate and interpolate clear sky fraction
  csf<-ukcpfutu_rsw/csr
  csf<-splina(.is(csf))
  csf[csf<0]<-0
  csf[csf>1]<-0
  # Calculate hourly clear sky fraction
  csr<-hourlyclm(year, month, day, lats)
  csr<-rast(dailycla(csr, dim(r)[[2]]+2)) # function also works for hourly
  ext(csr)<-ext(mnlon,mxlon,mnlat,mxlat)
  crs(csr)<-"EPSG:4326"
  csr<-project(csr,r)
  swdn<-.rast(csf,r)*csr
  if (silent == FALSE) utils::setTxtProgressBar(pb,5)
  # diffuse radiation
  difr<-mask(.rast(hourlydifa(.is(swdn),year,month,day,ll$lats,ll$lons),r),r)
  if (silent == FALSE) utils::setTxtProgressBar(pb,7)
  # LW down
  lwdn<-.rast(lwrada(.is(ukcpfutu_skyem),.is(temp)),r)
  if (silent == FALSE) utils::setTxtProgressBar(pb,8)
  # Wind speed
  # ~~ Calculate wind direction
  uu<-project(rast(ukcp_futu$uas),r)
  vv<-project(rast(ukcp_futu$vas),r)
  wd<-(atan2(.is(uu),.is(vv))*180/pi+180)%%360
  if (silent == FALSE) utils::setTxtProgressBar(pb,9)
  # ~~ Calculate u and v wind vectors
  uu<-splina(.is(ukcpfutu_ws)*cos(wd*pi/180))
  if (silent == FALSE) utils::setTxtProgressBar(pb,10)
  vv<-splina(.is(ukcpfutu_ws)*sin(wd*pi/180))
  if (silent == FALSE) utils::setTxtProgressBar(pb,11)
  ws<-.rast(sqrt(uu^2+vv^2),r)
  if (silent == FALSE) utils::setTxtProgressBar(pb,12)
  wd<-.rast((atan2(uu,vv)*180/pi+180)%%360,r)
  if (silent == FALSE) utils::setTxtProgressBar(pb,13)
  # precipitation
  prec<-.rast(precha(.is(ukcpfutu_prec)),r)
  if (silent == FALSE) utils::setTxtProgressBar(pb,14)
  # add times
  tmemn<-tmed[1]
  tmemx<-tmed[length(tmed)]+23*3600
  tme<-as.POSIXlt(seq(tmemn,tmemx,by=3600),tz="UTC")
  # create output list
  out<-list(temp=temp,relhum=rh,pres=pk,swdown=swdn,difrad=difr,lwdown=lwdn,windspeed=ws,winddir=wd,precip=prec)
  for (i in 1:9) {
    time(out[[i]])<-tme
    out[[i]]<-wrap(out[[i]])
  }
  if (silent == FALSE) utils::setTxtProgressBar(pb,15)
  return(out)
}
