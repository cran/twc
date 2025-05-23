#' Show data content
#'
#' The function \code{infoNC} displays the specification of the desired file.
#'
#' @importFrom raster brick
#' @importFrom methods setGeneric setMethod
#' @importFrom utils capture.output
#' @param x Raster* Object; character
#' @return character
#' @export

setGeneric("infoNC", function(x) standardGeneric("infoNC"))

#' @rdname infoNC
#' @method infoNC Raster

setMethod("infoNC", "Raster",
          function(x) {
            dummie <- capture.output(x)
            dummie <- dummie[-length(dummie)]
            return(dummie)
          })

#' @rdname infoNC
#' @method infoNC character

setMethod("infoNC", "character",
          function(x) {
            x <- brick(x)
            dummie <- capture.output(x)
            dummie <- dummie[-length(dummie)]
            return(dummie)
          })
