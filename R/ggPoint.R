##########################################################################################
# ggplot2 Wrapper Methods For Easy Plotting
##########################################################################################

#' A ggplot-based dot plot wrapper function
#'
#' This function is a wrapper around ggplot geom_point to allow for a more intuitive plotting of ArchR data.
#'
#' @param x A numeric vector containing the x-axis values for each point.
#' @param y A numeric vector containing the y-axis values for each point.
#' @param color A numeric/categorical vector used to determine the coloration for each point.
#' @param discrete A boolean value indicating whether the supplied data is discrete (`TRUE`) or continuous (`FALSE`).
#' @param discreteSet The name of a custom palette from `ArchRPalettes` to use for categorical/discrete color.
#' This argument is only used if `discrete` is set to `TRUE`.
#' @param continuousSet The name of a custom palette from `ArchRPalettes` to use for numeric color.
#' This argument is only used if `discrete` is set to `FALSE`.
#' @param labelMeans A boolean value indicating whether the mean of each categorical/discrete color should be labeled.
#' @param pal A custom palette used to override discreteSet/continuousSet for coloring vector.
#' @param defaultColor The default color for points that do not have another color applied (i.e. `NA` values).
#' @param highlightPoints A integer vector describing which points to hightlight. The remainder of points will be colored light gray.
#' @param colorDensity A boolean value indicating whether the density of points on the plot should be indicated by color.
#' If `TRUE`, continuousSet is used as the color palette.
#' @param size The numeric size of the points to be plotted.
#' @param xlim A numeric vector of two values indicating the lower and upper bounds of the x-axis on the plot.
#' @param ylim A numeric vector of two values indicating the lower and upper bounds of the y-axis on the plot.
#' @param extend A numeric value indicating the fraction to extend the x-axis and y-axis beyond the maximum and minimum
#' values if `xlim` and `ylim` are not provided. For example, 0.05 will extend the x-axis and y-axis by 5 percent on each end.
#' @param xlabel The label to plot for the x-axis.
#' @param ylabel The label to plot for the y-axis.
#' @param title The title of the plot.
#' @param randomize A boolean value indicating whether to randomize the order of the points when plotting.
#' @param seed A numeric seed number for use in randomization.
#' @param colorTitle A title to be added to the legend if `color` is supplied.
#' @param colorOrder A vector that allows you to control the order of palette colors associated with the values in `color`.
#' For example if you have `color` as `c("a","b","c")` and want to have the first color selected from the palette be used for
#' "c", the second color for "b", and the third color for "a", you would supply the `colorOrder` as `c("c", "b", "a")`.
#' @param colorLimits A numeric vector of two values indicating the lower and upper bounds of colors if numeric. Values
#' beyond these limits are thresholded.
#' @param alpha A number indicating the transparency to use for each point. See `ggplot2` for more details.
#' @param baseSize The base font size (in points) to use in the plot.
#' @param legendSize The size in inches to use for plotting the color legend.
#' @param ratioYX The aspect ratio of the x and y axes on the plot.
#' @param labelAsFactors A boolean indicating whether to label the `color` input as a numeric factor (`TRUE`) or with a character string (`FALSE`).
#' @param fgColor The foreground color of the plot.
#' @param bgColor The background color of the plot.
#' @param bgWidth The background relative width size of the halos in the labeling.
#' @param labelSize The numeric font size of labels.
#' @param addFit A string indicating the method to use for adding a fit/regression line to the plot (see `ggplot2::geom_smooth()` methods).
#' If set to `NULL`, no fit/regression line is added.
#' @param rastr A boolean value that indicates whether the plot should be rasterized using `ggrastr`. This does not rasterize
#' lines and labels, just the internal portions of the plot.
#' @param dpi The resolution in dots per inch to use for the plot.
#' @export
ggPoint <- function(
    x = NULL, 
    y = NULL, 
    color = NULL, 
    discrete = TRUE, 
    discreteSet = "stallion",
    continuousSet = "solarExtra", 
    labelMeans = TRUE,  
    pal = NULL, 
    defaultColor = "lightGrey",
    highlightPoints = NULL,
    colorDensity = FALSE,
    size = 1, 
    xlim = NULL, 
    ylim = NULL, 
    extend = 0.05, 
    xlabel = "x", 
    ylabel = "y", 
    title = "", 
    randomize = FALSE, 
    seed = 1,
    colorTitle = NULL, 
    colorOrder = NULL, 
    colorLimits = NULL,
    alpha = 1, 
    baseSize = 10, 
    legendSize = 3,
    ratioYX = 1, 
    labelAsFactors = TRUE,
    fgColor = "black", 
    bgColor = "white", 
    bgWidth = 1,
    labelSize = 3,
    addFit = NULL, 
    rastr = FALSE, 
    dpi = 300,
    ...
){

  stopifnot(length(y) == length(x))
  if(length(x) < 5){
    stop("x must be at least length 5 to plot!")
  }
  
  if(randomize){
    set.seed(seed)
    idx <- sample(seq_along(x), length(x))
  }else{
    idx <- seq_along(x)
  }
  
  df <- data.frame(x = x, y = y)
  include <- which(is.finite(x) & is.finite(y))
  
  if(length(include) != length(x)){
    message("Some values are not finite! Excluding these points!")
    df <- df[include,]
    x <- x[include]
    y <- y[include]
    if(!is.null(color)){
      color <- color[include]
    }
  }
  
  if(is.null(xlim)){
    xlim <- range(df$x) %>% extendrange(f = extend)
  }
  
  if(is.null(ylim)){
    ylim <- range(df$y) %>% extendrange(f = extend)
  }
  
  ratioXY <- ratioYX * diff(xlim)/diff(ylim)
  
  #Plot
  #.requirePackage("ggplot2", source = "cran")
  
  if (is.null(color) & !colorDensity) {
    
    p <- ggplot(df[idx,], aes(x = x, y = y)) + 
      coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = F) + 
      xlab(xlabel) + ylab(ylabel) + 
      ggtitle(title) + 
      theme_ArchR(baseSize = baseSize)
    
    if(rastr){
      p <- p + .geom_point_rast2(
        size = size, raster.dpi = dpi, alpha = alpha, color = defaultColor)
      # if(!requireNamespace("ggrastr", quietly = TRUE)){
      #   message("ggrastr is not available for rastr of points, continuing without rastr!")
      #   p <- p + geom_point(size = size, alpha = alpha, color = defaultColor)
      # }else{
      #   .requirePackage("ggrastr")
      #   p <- p + geom_point_rast(
      #       size = size, raster.dpi = dpi, alpha = alpha, color = defaultColor)
      # }
    }else{
      p <- p + geom_point(size = size, alpha = alpha, color = defaultColor)
    }
    
  }else {
    
    if(colorDensity){
      
      discrete <- FALSE
      df <- .getDensity(x, y, n = 100, sample = NULL) #change
      df <- df[order(df$density), ,drop=FALSE]
      df$color <- df$density
      
      if(is.null(colorTitle)){
        colorTitle <- "density"
      }
      
    }else if(discrete){
      
      if(!is.null(highlightPoints)){
        if(length(highlightPoints) < length(color)){
          color[-highlightPoints] <- "Non.Highlighted"
          idx <- c(idx[-highlightPoints], idx[highlightPoints])
        }
      }
      color <- paste0(color)
      
      if(!is.null(colorOrder)){
        if(!all(color %in% colorOrder)){
          stop("Not all colors are in colorOrder!")
        }
      }else{
        colorOrder <- gtools::mixedsort(unique(color))
      }
      
      if(is.null(colorTitle)){
        colorTitle <- "color"
      }
      
      stopifnot(length(color) == nrow(df))
      df$color <- factor(color, levels = colorOrder)
      
      if(labelAsFactors){
        df$color <- factor(
          x = paste0(paste0(match(paste0(df$color), paste0(levels(df$color)))), "-", paste0(df$color)), 
          levels = paste0(seq_along(levels(df$color)), "-", levels(df$color))
        )
        if(!is.null(pal)){
          #print(pal)
          #print(paste0(levels(df$color))[match(names(pal), colorOrder)])
          names(pal) <- paste0(levels(df$color))[match(names(pal), colorOrder)]
        }
        colorOrder <- paste0(levels(df$color))
      }
      
    }else{
      stopifnot(length(color) == nrow(df))
      if(!is.null(highlightPoints)){
        if(length(highlightPoints) < length(color)){
          color[-highlightPoints] <- NA
          idx <- c(idx[-highlightPoints], idx[highlightPoints])
        }
      }
      if(!is.null(colorLimits)){
        color[color < min(colorLimits)] <- min(colorLimits)
        color[color > max(colorLimits)] <- max(colorLimits)
      }
      df$color <- color
    }
    
    p <- ggplot(df[idx,], aes(x = x, y = y, color = color)) +  
      coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = FALSE) + 
      xlab(xlabel) + ylab(ylabel) + 
      ggtitle(title) + theme_ArchR(baseSize = baseSize) +
      theme(legend.direction = "horizontal", legend.box.background = element_rect(color = NA)) +
      labs(color = colorTitle)
    
    if(rastr){
      
      p <- p + .geom_point_rast2(
        size = size, raster.dpi = dpi, alpha = alpha, 
        raster.width = min(par('fin')), 
        raster.height = (ratioYX * min(par('fin')))
      )
      
      # if(!requireNamespace("ggrastr", quietly = TRUE)){
      #   message("ggrastr is not available for rastr of points, continuing without rastr!")
      #   message("To install ggrastr try : devtools::install_github('VPetukhov/ggrastr')")
      #   p <- p + geom_point(size = size, alpha = alpha)
      # }else{
      #   .requirePackage("ggrastr", installInfo = "devtools::install_github('VPetukhov/ggrastr')")
      #   p <- p + geom_point_rast(
      #       size = size, raster.dpi = dpi, alpha = alpha, 
      #       raster.width=par('fin')[1], 
      #       raster.height = (ratioYX * par('fin')[2])
      #     )
      # }
      
    }else{
      
      p <- p + geom_point(size = size, alpha = alpha)
      
    }
    
    if (discrete) {
      
      if (!is.null(pal)) {
        p <- p + scale_color_manual(values = pal)
      }else {
        pal <- paletteDiscrete(set = discreteSet, values = colorOrder)
        if(!is.null(highlightPoints)){
          pal[grep("Non.Highlighted", names(pal))] <- "lightgrey"
        }
        #print(pal)
        p <- p + scale_color_manual(values = pal) +
          guides(color = guide_legend(override.aes = list(size = legendSize, shape = 15)))
      }
      
      if (labelMeans) {
        
        dfMean <- split(df, df$color) %>% lapply(., function(x) {
          data.frame(x = median(x[, 1]), y = median(x[, 2]), color = x[1, 3])
        }) %>% Reduce("rbind", .)
        
        if(labelAsFactors){
          dfMean$label <- stringr::str_split(paste0(seq_len(nrow(dfMean))), pattern = "\\-", simplify=TRUE)[,1]
        }else{
          dfMean$label <- dfMean$color
        }
        dfMean$text <- stringr::str_split(dfMean$color, pattern = "-", simplify = TRUE)[,1]
        
        # make halo layers, similar to https://github.com/GuangchuangYu/shadowtext/blob/master/R/shadowtext-grob.R#L43
        theta <- seq(pi / 8, 2 * pi, length.out = 16)
        xo <- bgWidth * diff(range(df$x)) / 300
        yo <- bgWidth * diff(range(df$y)) / 300
        for (i in theta) {
          p <- p + 
            geom_text(data = dfMean, 
                      aes_q(
                        x = bquote(x + .(cos(i) * xo)),
                        y = bquote(y + .(sin(i) * yo)),
                        label = ~text
                      ),
                      size = labelSize,
                      color = bgColor
            )
        }
        
        if(is.null(fgColor)){
          p <- p + geom_text(data = dfMean, aes(x = x, y = y, color = color, label = label), size = labelSize, show.legend = FALSE)
        }else{
          p <- p + geom_text(data = dfMean, aes(x = x, y = y, label = label), color = fgColor, size = labelSize, show.legend = FALSE) 
        }
        
      }
      
    }else{
      
      if (!is.null(pal)) {
        if(!is.null(colorLimits)){
          p <- p + scale_colour_gradientn(colors = pal, limits=colorLimits, na.value = "lightgrey")
        }else{
          p <- p + scale_colour_gradientn(colors = pal, na.value = "lightgrey")
        }
      }else {
        if(!is.null(colorLimits)){
          p <- p + scale_colour_gradientn(colors = paletteContinuous(set = continuousSet), limits=colorLimits, na.value = "lightgrey")
        }else{
          p <- p + scale_colour_gradientn(colors = paletteContinuous(set = continuousSet), na.value = "lightgrey")
        }
      }
    }
    
  }
  
  if (!is.null(addFit)) {
    p <- p + geom_smooth(data = df, aes(color = NULL), method = addFit, color = "black") + 
      ggtitle(paste0(title, "\nPearson = ", round(cor(df$x, df$y), 3), "\nSpearman = ", round(cor(df$x, df$y, method = "spearman"), 3)))
  }
  
  p <- p + theme(legend.position = "bottom", legend.key = element_rect(size = 2))#, legend.spacing.x = unit(0.1, 'cm'), legend.spacing.y = unit(0.1, 'cm'))
  
  if(!is.null(ratioYX)){
    attr(p, "ratioYX") <- ratioYX
  }
  
  return(p)
  
}



#' ggplot2 default theme for ArchR
#'
#' This function returns a ggplot2 theme that is black borded with black font.
#' 
#' @param color The color to be used for text, lines, ticks, etc for the plot.
#' @param textFamily The font default family to be used for the plot.
#' @param baseSize The base font size (in points) to use in the plot.
#' @param baseLineSize The base line width (in points) to be used throughout the plot.
#' @param baseRectSize The base line width (in points) to use for rectangular boxes throughout the plot.
#' @param plotMarginCm The width in centimeters of the whitespace margin around the plot.
#' @param legendPosition The location to put the legend. Valid options are "bottom", "top", "left", and "right.
#' @param legendTextSize The base text size (in points) for the legend text.
#' @param axisTickCm The length in centimeters to be used for the axis ticks.
#' @param xText90 A boolean value indicating whether the x-axis text should be rotated 90 degrees counterclockwise.
#' @param yText90 A boolean value indicating whether the y-axis text should be rotated 90 degrees counterclockwise.
#' @export
theme_ArchR <- function(
    color = "black",
    textFamily = "sans",
    baseSize = 10, 
    baseLineSize = 0.5,
    baseRectSize = 0.5,
    plotMarginCm = 1,
    legendPosition = "bottom",
    legendTextSize = 5,
    axisTickCm = 0.1,
    xText90 = FALSE,
    yText90 = FALSE
){

  theme <- theme_bw() + theme(
    text = element_text(family = textFamily),
    axis.text = element_text(color = color, size = baseSize), 
    axis.title = element_text(color = color, size = baseSize),
    title = element_text(color = color, size = baseSize),
    plot.margin = unit(c(plotMarginCm, plotMarginCm, plotMarginCm, plotMarginCm), "cm"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = color, size = (4/3) * baseRectSize * as.numeric(grid::convertX(grid::unit(1, "points"), "mm"))),
    axis.ticks.length = unit(axisTickCm, "cm"), 
    axis.ticks = element_line(color = color, size = baseLineSize * (4/3) * as.numeric(grid::convertX(grid::unit(1, "points"), "mm"))),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.text = element_text(color = color, size = legendTextSize),
    legend.box.background = element_rect(color = NA),
    #legend.box.background = element_rect(fill = "transparent"),
    legend.position = legendPosition,
    strip.text = element_text(size = baseSize, color="black")#,
    #plot.background = element_rect(fill = "transparent", color = NA)
  )
  
  if(xText90){
    theme <- theme %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  
  if(yText90){
    theme <- theme %+replace% theme(axis.text.y = element_text(angle = 90, vjust = 1))
  }
  
  return(theme)
  
}



##########################################################################################
# ggplot2 helper functions
##########################################################################################

.checkCairo <- function(){
  tryCatch({
    tmp <- dev.cur()
    Cairo::Cairo(type='raster')
    dev.off()
    dev.set(tmp)
    TRUE
  }, error = function(e){
    FALSE
  })
}

## Adapted from 
## https://github.com/tidyverse/ggplot2/blob/660aad2db2b3495ae0d8040915a40d247133ffc0/R/geom-point.r
## from https://github.com/VPetukhov/ggrastr/blob/master/R/geom-point-rast.R
## This funciton now handles issues with Cairo installation that can lead to plot errors
.geom_point_rast2 <- function(
    mapping = NULL,
    data = NULL,
    stat = "identity",
    position = "identity",
    ...,
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE,
    raster.width = min(par('fin')), 
    raster.height = min(par('fin')), 
    raster.dpi = 300
){
  
  GeomPointRast <- tryCatch({
    
    if(!.checkCairo()){
      stop()
    }
    
    #Try to create a geom rast for points if not then just use normal geom_point
    ggplot2::ggproto(
      "GeomPointRast",
      ggplot2::GeomPoint,
      required_aes = c("x", "y"),
      non_missing_aes = c("size", "shape", "colour"),
      default_aes = aes(
        shape = 19, colour = "black", size = 1.5, fill = NA,
        alpha = NA, stroke = 0.5
      ),
      
      draw_panel = function(data, panel_params, coord, na.rm = FALSE, 
                            raster.width=min(par('fin')), raster.height=min(par('fin')), raster.dpi=300){
        
        #From ggrastr  
        prevDevID <- dev.cur()
        
        p <- ggplot2::GeomPoint$draw_panel(data, panel_params, coord)
        
        devID <- Cairo::Cairo(
          type='raster', 
          width=raster.width*raster.dpi, 
          height=raster.height*raster.dpi, 
          dpi=raster.dpi, 
          units='px', 
          bg="transparent"
        )[1]
        
        grid::pushViewport(grid::viewport(width=1, height=1))
        
        grid::grid.points(
          x=p$x, 
          y=p$y, 
          pch = p$pch, 
          size = p$size,
          name = p$name, 
          gp = p$gp, 
          vp = p$vp, 
          draw = TRUE
        )
        
        grid::popViewport()
        gridCapture <- grid::grid.cap()
        
        dev.off(devID)
        
        dev.set(prevDevID)
        
        grid::rasterGrob(
          gridCapture, 
          x=0, 
          y=0, 
          width = 1,
          height = 1,
          default.units = "native",
          just = c("left","bottom")
        )
        
      }
      
    )
    
  }, error = function(e){
    
    if(.checkCairo()){
      message("WARNING: Error found with trying to rasterize geom. Continuing without rasterization.")
    }else{
      message("WARNING: Error found with Cairo installation. Continuing without rasterization.")
    }
    
    #Default geom_point
    ggplot2::ggproto(
      "GeomPoint", 
      ggplot2::GeomPoint,
      required_aes = c("x", "y"),
      non_missing_aes = c("size", "shape", "colour"),
      default_aes = aes(
        shape = 19, colour = "black", size = 1.5, fill = NA,
        alpha = NA, stroke = 0.5
      ),
      
      draw_panel = function(data, panel_params, coord, na.rm = FALSE, 
                            raster.width=min(par('fin')), raster.height=min(par('fin')), raster.dpi=300){
        if (is.character(data$shape)) {
          data$shape <- ggplot2:::translate_shape_string(data$shape) #Hidden ggplot2
        }
        
        coords <- coord$transform(data, panel_params)
        
        pGrob <- grid::pointsGrob(
          x = coords$x, 
          y = coords$y,
          pch = coords$shape,
          gp = grid::gpar(
            col = scales::alpha(coords$colour, coords$alpha),
            fill = scales::alpha(coords$fill, coords$alpha),
            # Stroke is added around the outside of the point
            fontsize = coords$size * .pt + coords$stroke * .stroke / 2,
            lwd = coords$stroke * .stroke / 2
          )
        )
        
        pGrob
        
      },
      
      draw_key = ggplot2::draw_key_point
    )
    
    
  })
  
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomPointRast,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      raster.width=raster.width,
      raster.height=raster.height,
      raster.dpi=raster.dpi,
      ...
    )
  )
  
}



##########################################################################################
# Plot Aesthetics Objects and Methods
##########################################################################################

#' List of color palettes that can be used in plots
#' 
#' A collection of some original and some borrowed color palettes to provide appealing color aesthetics for plots in ArchR
#' 
#' @export
ArchRPalettes <- list(
  
  #DISCLOSURE: This is a collection of palettes that includes some original palettes and some palettes originally
  #implemented by others in other packages.
  #They are included here for convenience because they help improve plot aesthetics.
  
  #NOTE: all palettes included in the "Primarily Continuous Palettes" section should also work for discrete usage but not vice versa.
  #Each continuous palette has been ordered by color to generate a visually appealing discrete palette.
  
  #---------------------------------------------------------------
  # Primarily Discrete Palettes
  #---------------------------------------------------------------
  
  #20-colors
  stallion = c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
               "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D"),
  
  stallion2 = c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
                "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767"),
  
  calm = c("1"="#7DD06F", "2"="#844081", "3"="#688EC1", "4"="#C17E73", "5"="#484125", "6"="#6CD3A7", "7"="#597873","8"="#7B6FD0", "9"="#CF4A31", "10"="#D0CD47",
           "11"="#722A2D", "12"="#CBC594", "13"="#D19EC4", "14"="#5A7E36", "15"="#D4477D", "16"="#403552", "17"="#76D73C", "18"="#96CED5", "19"="#CE54D1", "20"="#C48736"),
  
  kelly = c("1"="#FFB300", "2"="#803E75", "3"="#FF6800", "4"="#A6BDD7", "5"="#C10020", "6"="#CEA262", "7"="#817066", "8"="#007D34", "9"="#F6768E", "10"="#00538A",
            "11"="#FF7A5C", "12"="#53377A", "13"="#FF8E00", "14"="#B32851", "15"="#F4C800", "16"="#7F180D", "17"="#93AA00", "18"="#593315", "19"="#F13A13", "20"="#232C16"),
  
  #16-colors
  bear = c("1"="#faa818", "2"="#41a30d","3"="#fbdf72", "4"="#367d7d",  "5"="#d33502", "6"="#6ebcbc", "7"="#37526d",
           "8"="#916848", "9"="#f5b390", "10"="#342739", "11"="#bed678","12"="#a6d9ee", "13"="#0d74b6",
           "14"="#60824f","15"="#725ca5", "16"="#e0598b"),
  
  #15-colors
  ironMan = c("9"='#371377',"3"='#7700FF',"2"='#9E0142',"10"='#FF0080', "14"='#DC494C',"12"="#F88D51","1"="#FAD510","8"="#FFFF5F","4"='#88CFA4',
              "13"='#238B45',"5"="#02401B", "7"="#0AD7D3","11"="#046C9A", "6"="#A2A475", "15"='grey35'),
  
  circus = c("1"="#D52126", "2"="#88CCEE", "3"="#FEE52C", "4"="#117733", "5"="#CC61B0", "6"="#99C945", "7"="#2F8AC4", "8"="#332288",
             "9"="#E68316", "10"="#661101", "11"="#F97B72", "12"="#DDCC77", "13"="#11A579", "14"="#89288F", "15"="#E73F74"),
  
  #12-colors
  paired = c("9"="#A6CDE2","1"="#1E78B4","3"="#74C476","12"="#34A047","11"="#F59899","2"="#E11E26",
             "10"="#FCBF6E","4"="#F47E1F","5"="#CAB2D6","8"="#6A3E98","6"="#FAF39B","7"="#B15928"),
  
  #11-colors
  grove = c("11"="#1a1334","9"="#01545a","1"="#017351","6"="#03c383","8"="#aad962","2"="#fbbf45","10"="#ef6a32","3"="#ed0345","7"="#a12a5e","5"="#710162","4"="#3B9AB2"),
  
  #7-colors
  summerNight = c("1"="#2a7185", "2"="#a64027", "3"="#fbdf72","4"="#60824f","5"="#9cdff0","6"="#022336","7"="#725ca5"),
  
  #5-colors
  zissou = c("1"="#3B9AB2", "4"="#78B7C5", "3"="#EBCC2A", "5"="#E1AF00", "2"="#F21A00"), #wesanderson
  darjeeling = c("1"="#FF0000", "2"="#00A08A", "3"="#F2AD00", "4"="#F98400", "5"="#5BBCD6"), #wesanderson
  rushmore = c("1"="#E1BD6D", "5"="#EABE94", "2"="#0B775E", "4"="#35274A" , "3"="#F2300F"), #wesanderson
  captain = c("1"="grey","2"="#A1CDE1","3"="#12477C","4"="#EC9274","5"="#67001E"),
  
  #---------------------------------------------------------------
  # Primarily Continuous Palettes
  #---------------------------------------------------------------
  
  #10-colors
  horizon = c("1"='#000075',"4"='#2E00FF', "6"='#9408F7', "10"='#C729D6', "8"='#FA4AB5', "3"='#FF6A95', "7"='#FF8B74', "5"='#FFAC53', "9"='#FFCD32', "2"='#FFFF60'),
  
  #9-colors
  horizonExtra =c("1"="#000436","4"="#021EA9","6"="#1632FB","8"="#6E34FC","3"="#C732D5","9"="#FD619D","7"="#FF9965","5"="#FFD32B","2"="#FFFC5A"),
  blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D"),
  sambaNight = c("6"='#1873CC',"2"='#1798E5',"8"='#00BFFF',"5"='#4AC596',"1"='#00CC00',"4"='#A2E700',"9"='#FFFF00',"7"='#FFD200',"3"='#FFA500'), #buencolors
  solarExtra = c("5"='#3361A5', "7"='#248AF3', "1"='#14B3FF', "8"='#88CEEF', "9"='#C1D5DC', "4"='#EAD397', "3"='#FDB31A',"2"= '#E42A2A', "6"='#A31D1D'),  #buencolors
  whitePurple = c("9"='#f7fcfd',"6"='#e0ecf4',"8"='#bfd3e6',"5"='#9ebcda',"2"='#8c96c6',"4"='#8c6bb1',"7"='#88419d',"3"='#810f7c',"1"='#4d004b'),
  whiteBlue = c("9"='#fff7fb',"6"='#ece7f2',"8"='#d0d1e6',"5"='#a6bddb',"2"='#74a9cf',"4"='#3690c0',"7"='#0570b0',"3"='#045a8d',"1"='#023858'),
  whiteRed = c("1"="white", "2"="red"),
  comet = c("1"="#E6E7E8","2"="#3A97FF","3"="#8816A7","4"="black"),
  
  #7-colors
  greenBlue = c("4"='#e0f3db',"7"='#ccebc5',"2"='#a8ddb5',"5"='#4eb3d3',"3"='#2b8cbe',"6"='#0868ac',"1"='#084081'),
  
  #6-colors
  beach = c("4"="#87D2DB","1"="#5BB1CB","6"="#4F66AF","3"="#F15F30","5"="#F7962E","2"="#FCEE2B"),
  
  #5-colors
  coolwarm = c("1"="#4858A7", "4"="#788FC8", "5"="#D6DAE1", "3"="#F49B7C", "2"="#B51F29"),
  fireworks = c("5"="white","2"="#2488F0","4"="#7F3F98","3"="#E22929","1"="#FCB31A"),
  greyMagma = c("2"="grey", "4"="#FB8861FF", "5"="#B63679FF", "3"="#51127CFF", "1"="#000004FF"),
  fireworks2 = c("5"="black", "2"="#2488F0","4"="#7F3F98","3"="#E22929","1"="#FCB31A"),
  purpleOrange = c("5"="#581845", "2"="#900C3F", "4"="#C70039", "3"="#FF5744", "1"="#FFC30F")
)

#' Optimized discrete color palette generation
#'
#' This function assesses the number of inputs and returns a discrete color palette that is tailored to provide the most
#' possible color contrast from the designated color set.
#'
#' @param values A character vector containing the sample names that will be used. Each entry in this character vector will be
#' given a unique color from the designated palette set.
#' @param set The name of a color palette provided in the `ArchRPalettes` list object.
#' @param reverse A boolean variable that indicates whether to return the palette colors in reverse order.
#' @export
paletteDiscrete <- function(
    values = NULL,
    set = "stallion",  
    reverse = FALSE
){
  

  values <- unique(values)
  values <- gtools::mixedsort(values)
  n <- length(unique(values))
  pal <- ArchRPalettes[[set]]
  palOrdered <- pal[gtools::mixedsort(names(pal))] #mixed sort gets 1,2,3,4..10,11,12
  
  if(n > length(palOrdered)){
    message("Length of unique values greater than palette, interpolating..")
    palOut <- colorRampPalette(pal)(n)
  }else{
    palOut <- palOrdered[seq_len(n)]
  }
  
  if(reverse){
    palOut <- rev(palOut)
  }
  
  names(palOut) <- unique(values)
  
  return(palOut)
  
}

#' Continuous Color Palette
#'
#' @param set The name of a color palette provided in the `ArchRPalettes` list object.
#' @param n The number of unique colors to generate as part of this continuous color palette.
#' @param reverse A boolean variable that indicates whether to return the palette colors in reverse order.
#' @export
paletteContinuous <- function(
    set = "solarExtra", 
    n = 256, 
    reverse = FALSE
){
  

  pal <- ArchRPalettes[[set]]
  palOut <- colorRampPalette(pal)(n)
  
  if(reverse){
    palOut <- rev(palOut)
  }
  
  return(palOut)
  
}