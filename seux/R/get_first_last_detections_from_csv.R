#' @title Get first and last detections from csv file
#'
#' @description Returns the first and last detection years
#'  stored in a csv file after performing basic checks.
#'
#' @param fName String, a filename.
#'
#' @param frst_col Integer or String, the column containing each species' year of first detection.
#'
#' @param last_col Integer or String, the column containing each species' year of last detection.
#'
#' @return A dataframe with two columns: frstDetn, lastDetn
#'  \describe{
#'   \item{frstDetn:}{Integers, species' year of first detection.}
#'   \item{lastDetn:}{Integers, species' year of last detection.}
#'  }
#'
#' @importFrom utils read.csv
#'
#' @export

get_first_last_detections_from_csv <- function(fName, frst_col, last_col) {

    # check file name

    if ( ! is.character(fName) ) { stop("Input fName needs to be a string giving the file name: ", as.character(fName)) }

    # read in data

    frst_last <- read.csv(fName,header=T)

    # check headers

    if ( is.character(frst_col) ) {
        frst_col <- gsub(" ", ".", frst_col)
        if (! frst_col %in% colnames(frst_last) ) { 
            stop( "Could not find header in csv file: ", as.character(frst_col) )
        } else {
            frstDetn <- frst_last[frst_col][,1]
        }
    } else if ( is.numeric(frst_col) ) {
        if ( frst_col %% 1 == 0 ) { # if it's an whole number
            frstDetn <- frst_last[,frst_col] 
        } else {
            stop( "Column for first-detection dates needs to be whole no. or string of header: ", as.character(frst_col) )
        }
    } else {
        stop( "Column for first-detection dates needs to be whole no. or string of header." )
    }

    if ( is.character(last_col) ) {
        last_col <- gsub(" ", ".", last_col)
        if (! last_col %in% colnames(frst_last) ) { 
            stop( "Could not find header in csv file: ", as.character(last_col) )
        } else {
            lastDetn <- frst_last[last_col][,1]
        }
    } else if ( is.numeric(last_col) ) {
        if ( last_col %% 1 == 0 ) { # if it's an whole number
            lastDetn <- frst_last[,last_col] 
        } else {
            stop( "Column for last-detection dates needs to be whole no. or string of header: ", as.character(last_col) )
        }
    } else {
        stop( "Column for last-detection dates needs to be whole no. or string of header." )
    }

    # check data

    input_okay <- check_frst_last_detns(frstDetn, lastDetn)
    if (! input_okay$flag ){
        stop( input_okay$message )
    }

    return( data.frame( frstDetn=frstDetn, lastDetn=lastDetn ) )

}
