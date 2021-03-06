#' plot clinical information
#'
#' given a data frame with columns names sample, variable, and value create a
#' ggplot2 object
#' @name multi_buildClin
#' @param x a data frame in "long" format giving additional information to be
#' plotted, requires columns "sample", "variable", and "value"
#' @param clin.legend.col an integer specifying the number of columns to plot in
#' the legend
#' @param clin.var.colour a named character vector specifying the mapping
#' between colors and variables
#' @param clin.var.order a character vector of variables to order the legend by
#' @param clin.layers additional ggplot2 layers to plot
#' @return a grob object
#' @import ggplot2
#' @importFrom viridis viridis
#' @noRd

multi_buildClin <- function(x, clin.legend.col=1, clin.var.colour=NULL,
                                clin.var.order=NULL, clin.layers=NULL)
{
    # Define parameters
    x_label <- xlab(paste0("Sample n=", length(unique(x$sample))))
    leg_guide <- guides(fill=guide_legend(ncol=clin.legend.col))
    
    # Set the clinical variable order and colour
    if(!is.null(clin.var.colour))
    {
        na_index <- which(grepl("^NA$", names(clin.var.colour)))
        # check that there is only na value
        if(length(na_index) == 1)
        {
            na_colour <- clin.var.colour[na_index]
            clin.var.colour <- clin.var.colour[-na_index]
        } else {
            na_colour <- NA
        }
    } else {
        na_colour <- NA
    }

    if(!is.null(clin.var.colour) & is.null(clin.var.order))
    {
        clin_fill_colour <- scale_fill_manual(name="Clinical Data",
                                              values=clin.var.colour,
                                              na.value=na_colour)
    } else if(!is.null(clin.var.colour) & !is.null(clin.var.order)) {
        clin_fill_colour <- scale_fill_manual(name="Clinical Data",
                                              breaks=clin.var.order,
                                              values=clin.var.colour,
                                              na.value=na_colour)
    } else if(is.null(clin.var.colour) & !is.null(clin.var.order)) {
        clin_fill_colour <- scale_fill_manual(name="Clinical Data",
                                              breaks=clin.var.order,
                                              values=viridis::viridis(length(unique(x$value))),
                                              na.value=na_colour)
    } else {
        clin_fill_colour <- geom_blank()
    }

    if(!is.null(clin.layers))
    {
        layers <- clin.layers
    } else {
        layers <- geom_blank()
    }

    # Define the theme
    theme <- theme(panel.grid.major = element_blank(),
                   panel.grid.minor=element_blank(),
                   panel.background=element_rect(fill='white',
                                                 colour='white'),
                   axis.ticks.x=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.text.x=element_blank(),
                   axis.title.x=element_text(size=16),
                   legend.title=element_text(size=14),
                   axis.title.y=element_blank(),
                   axis.text.y=element_text(size=14, colour='black'),
                   legend.position='right')

    # Define the main plot
    p1 <- ggplot(x, aes_string(x='sample', y='variable', fill='value')) +
    geom_tile() + theme + x_label + leg_guide + clin_fill_colour + layers

    return(p1)
}
