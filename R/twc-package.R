#' twc: Terrestrial Water Cycle
#'
#' twc is an R package that that constitutes the core functions to analyze terrestrial water cycle data across various spatio-temporal scales.
#'
#' @author \strong{Coded by:} Mijael Rodrigo Vargas Godoy \email{vargas_godoy@fzp.czu.cz}
#' @author \strong{Conceptual design by:} Mijael Rodrigo Vargas Godoy \email{vargas_godoy@fzp.czu.cz} and Yannis Markonis \email{markonis@fzp.czu.cz}
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import data.table doParallel foreach ncdf4 parallel sp
#' @importFrom lubridate days_in_month
#' @importFrom magrittr %>%
#' @importFrom methods as is setGeneric setMethod
#' @importFrom raster aggregate area as.data.frame brick calc cellStats crop extent getValues getZ mask nlayers raster rasterFromXYZ res setZ subset xFromCol yFromRow zApply
#' @importFrom sf read_sf st_bbox
#' @importFrom utils capture.output tail
## usethis namespace: end
NULL
