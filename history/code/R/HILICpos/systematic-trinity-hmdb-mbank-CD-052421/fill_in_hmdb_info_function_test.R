#' This part retrieve metabolite information from hmdb based on HMDB0000001

.hmxPath = function(prefix="http://www.hmdb.ca/metabolites/",
                    id="HMDB0000001", ...) {
  sub("__PRE__", prefix, sub("%%ID%%", id, "__PRE__%%ID%%.xml") )
}

hmxToList = function (prefix = "http://www.hmdb.ca/metabolites/", id = "HMDB0000001", 
                      ...) 
{
  requireNamespace("XML")
  stopifnot(is.atomic(prefix), length(prefix)==1, is.atomic(id), length(id)==1)
  txt = readLines(.hmxPath(prefix=prefix, id=id, ...))
  prs = xmlTreeParse(txt, asText=TRUE)
  xmlToList(prs)
}

library("XML")
hmxToList(prefix="http://www.hmdb.ca/metabolites/",
          id="HMDB0014328")


#' @
#' @
#' @
#' @
#' @
  
#' This part attempts retrieve spectra information from hmdb based on spectra ID e.g., 2619.
##' However the output seems not be very helpful. 
##' And as I realized that in the below function, Johannes already found no precursorMz is given
###' https://github.com/EuracBiomedicalResearch/CompoundDb/blob/master/R/spectrum-import-functions.R
###' Thus I am not going to use this to calculate precursorMz. 
###' I will just access the HMDB and then assume that all positive mode will be M+H; negative M-H
###' And thus, I can calculate by adding on exact mass with plus/minus H.

# source: https://hmdb.ca/spectra/ms_ms/2619/generate_mzml

.hmxPath = function(prefix="http://www.hmdb.ca/spectra/ms_ms/",
                    id="", ...) {
  sub("__PRE__", prefix, sub("%%ID%%", id, "__PRE__%%ID%%/generate_mzml") )
}

hmxToList = function (prefix = "http://www.hmdb.ca/spectra/ms_ms/", id = "HMDB0000001", 
                      ...) 
{
  requireNamespace("XML")
  stopifnot(is.atomic(prefix), length(prefix)==1, is.atomic(id), length(id)==1)
  txt = readLines(.hmxPath(prefix=prefix, id=id, ...))
  prs = xmlTreeParse(txt, asText=TRUE)
  xmlToList(prs)
}

library("XML")
hmxToList(prefix="http://www.hmdb.ca/spectra/ms_ms/",
          id="248219")
