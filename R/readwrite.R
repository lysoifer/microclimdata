#' @title Write timeseries profiles from RunModelFull
#' @description
#' Write output from RunModelFull
#' @param mout output of RunModelFull
#' @param fname file name (must be .h5)
#' @param h heights of canopy nodes in mout (represent top of each bin)
#' @param st_time start time
#' @param en_time end time - assumes time step between start and end times is 1 hour
#' @param microout vector of parameters to save; can be any of the names of mout (e.g., tair, tleaf, relhum)
#' @param vegp if not NULL, writes vegp inputs
#' @param soilp if not NULL, writes soilc inputs
#' @param macroclim if not null writes macroclimate data inputs
#' @param zref if not NULL writes zref parameter
#' @param Lfrac if not NULL writes Lfrac inputs
#' @param paii if not NULL writes paii inputs
#' 
#' @returns writes an h5 file with three groups. 1) microclim/ contains modelled
#'  microclimate arrays. Writes the variables specified in microout. Columns represent
#'  time (hour), and rows represent height in the canopy. Columns are ordered from
#'  first(left) to last (right), and heights are ordered from top of canopy (1)
#'  to bottom of canopy (nrow (microclim)), Therefore microclim[[1]][1,1]
#'  represents the top of the canopy at the first time step and microclim[[1]][nrow(out),1]
#'  represents the bottom of the canopy at the first time step.
#'  The microclim group contains four attributes 1) height of each canopy node (top to bottom),
#'  start time (st_time), end time (en_time), and time step in hours. Time variables
#'  allow recreation of the obs_time vector used in the models.
#'  Group 2) inputs - if inputs are provided, they are written the the file
#'  Group 3) metadata - contains, lat, lon, time_zone
#'  @export
write_micropoint_full <- function(mout, fname, h, st_time, en_time, microout = c("tair", "tleaf", "relhum"), 
                                vegp = NULL, soilp = NULL, macroclim = NULL, 
                                zref = NULL, Lfrac = NULL, paii, overwrite = T) {
  
  
  h <- h[length(h):1] # flip h so top of matrix[1,1] is top of canopy
  
  if(overwrite) unlink(fname)
  
  # create file
  h5createFile(fname)
  
  # create a group
  h5createGroup(fname, "microclim")
  
  
  if("tair" %in% microout) {
    tair = mout[["tair"]]
    tair = t(tair)
    tair = tair[nrow(tair):1,]
    h5write(tair, file = fname, name = "microclim/tair")
  }
  
  if("tleaf" %in% microout) {
    tleaf = mout[["tleaf"]]
    tleaf = t(tleaf) # indexed from bottom (1) to top of canopy
    tleaf = tleaf[nrow(tleaf):1,] # now indexed from top (1) to bottom of canopy
    h5write(tleaf, file = fname, name = "microclim/tleaf")
  }
  
  if("relhum" %in% microout) {
    rh = mout[["relhum"]]
    rh = t(rh)
    rh = rh[nrow(rh):1,]
    h5write(rh, file = fname, name = "microclim/relhum")
  }
  
  if("Rdirdown" %in% microout) {
    Rdirdown = mout[["Rdirdown"]]
    Rdirdown = t(Rdirdown)
    Rdirdown = Rdirdown[nrow(Rdirdown):1,]
    h5write(Rdirdown, file = fname, name = "microclim/Rdirdown")
  }
  
  if("Rdifdown" %in% microout) {
    Rdifdown = mout[["Rdifdown"]]
    Rdifdown = t(Rdifdown)
    Rdifdown = Rdifdown[nrow(Rdifdown):1,]
    h5write(Rdifdown, file = fname, name = "microclim/Rdifdown")
  }
  
  if("Rswup" %in% microout) {
    Rswup = mout[["Rswup"]]
    Rswup = t(Rswup)
    Rswup = Rswup[nrow(Rswup):1,]
    h5write(Rswup, file = fname, name = "microclim/Rswup")
  }
  
  if("Rlwdown" %in% microout) {
    Rlwdown = mout[["Rlwdown"]]
    Rlwdown = t(Rlwdown)
    Rlwdown = Rlwdown[nrow(Rlwdown):1,]
    h5write(Rlwdown, file = fname, name = "microclim/Rlwdown")
  }
  
  if("Rlwup" %in% microout) {
    Rlwup = mout[["Rlwup"]]
    Rlwup = t(Rlwup)
    Rlwup = Rlwup[nrow(Rlwup):1,]
    h5write(Rlwup, file = fname, name = "microclim/Rlwup")
  }
  
  if("windspeed" %in% microout) {
    windspeed = mout[["windspeed"]]
    windspeed = t(windspeed)
    windspeed = windspeed[nrow(windspeed):1,]
    h5write(windspeed, file = fname, name = "microclim/windspeed")
  }
  
  if("tsoil" %in% microout) {
    tsoil = mout[["tsoil"]]
    tsoil = t(tsoil)
    tsoil = tsoil[nrow(tsoil):1,]
    h5write(tsoil, file = fname, name = "microclim/tsoil")
  }
  
  if("theta" %in% microout) {
    theta = mout[["theta"]]
    theta = t(theta)
    theta = theta[nrow(theta):1,]
    h5write(theta, file = fname, name = "microclim/theta")
  }
  
  h5createGroup(fname, "inputs/")
  if(!is.null(vegp)) {
    h5write(vegp, file = fname, name = "inputs/vegp")
  }
  if(!is.null(soilc)) {
    h5write(soilc, file = fname, name = "inputs/soilc")
  }
  if(!is.null(macroclim)) {
    h5write(macroclim, file = fname, name = "inputs/macroclim")
  }
  if(!is.null(zref)) {
    h5write(zref, file = fname, name = "inputs/zref")
  }
  if(!is.null(Lfrac)) {
    h5write(lfrac, file = fname, name = "inputs/Lfrac")
  }
  
  H5close()
  
  fid <- H5Fopen(fname)
  g <- H5Gopen(fid, "microclim/")
  h5writeAttribute(h, g, name = "height")
  h5writeAttribute(as.character(st_time), g, name = "start_time")
  h5writeAttribute(as.character(en_time), g, name = "end_time")
  h5writeAttribute(3600, g, name = "time_step_seconds")
  H5close()
  
  h5createGroup(fname, "metadata/")
  h5write(x, file = fname, name = "metadata/lon")
  h5write(y, file = fname, name = "metadata/lat")
  h5write("UTC", file = fname, name = "metadata/time_zone")
  h5closeAll()
} 
