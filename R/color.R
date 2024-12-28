#' Generate colors from a customed color palette. This function is from R package: CellChat.
#' 
#' @param n number of colors. Max number is 26.
#' 
#' @return A color palette for plotting
#' 
#' @export
scPalette1 <- function (n) 
{
  colorSpace <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                  "#F29403", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4", 
                  "#222F75", "#1B9E77", "#B2DF8A", "#E3BE00", "#FB9A99", 
                  "#E7298A", "#910241", "#00CDD1", "#A6CEE3", "#CE1261", 
                  "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", "#B3DE69", 
                  "#8DD3C7", "#999999")
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  }
  else {
    colors <- (grDevices::colorRampPalette(colorSpace))(n)
  }
  return(colors)
}

#' Generate colors from a customed color palette
#' 
#' @param n number of colors. Max number is 30.
#' 
#' @return A color palette for plotting
#' 
#' @export
scPalette2 <- function (n) 
{
  colorSpace <- c("#00A087FF", "#4DBBD5FF", "#E64B35FF", "#3C5488FF", 
                  "#F38400", "#A1CAF1", "#BE0032", "#C2B280", "pink", "#008856", 
                  "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C", 
                  "#DCD300", "#882D17", "#8DB600", "#654522", "#E25822", 
                  "#2B3D26", "#848482", "#008856", "#E68FAC", "#0067A5", 
                  "#604E97", "#F6A600", "#B3446C", "#DCD300", "#882D17")
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  }
  else {
    colors <- (grDevices::colorRampPalette(colorSpace))(n)
  }
  return(colors)
}

#' Generate ggplot2's default colors 
#' 
#' @param n number of colors. 
#' 
#' @return A color palette for plotting
#' 
#' @export
scPalette4 <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#' Generate colors from https://doi.org/10.1016/j.cell.2020.03.048, Fig1A.
#' 
#' @param n number of colors. Max number is 44.
#' 
#' @return A color palette for plotting
#' 
#' @export
scPalette3 <- function (n) 
{
  colorSpace <- c("#b589bc", "#ecb888", "#fc7440", "#a032cb", 
                  "#cd8ab1", "#e9d6eb", "#77a1ec", "#a3c4a5", "#eed785", 
                  "#d4b595", "#efe0e7", "#eebed6", "#f7d39b", "#cd7560", 
                  "#7fc28e", "#f9cabe", "#c76da8", "#fef6ae", "#6276b2", 
                  "#e69d8e", "#da7574", "#c5d9df", "#e16db7", "#908ebc", 
                  "#af88bb", "#dedbee", "#f07590", "#dfb9d5", "#b099b5", 
                  "#5394c3", "#b8a89f", "#917393", "#85c1d7", "#dbebf7", 
                  "#84d7f7", "#d4eef6", "#f5bcc9", "#5eb4f3", "#f7f4b7", 
                  "#6f89d3", "#50c2fb", "#ea94cd", "#d8b6f7", "#efc1ec")
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  }
  else {
    colors <- (grDevices::colorRampPalette(colorSpace))(n)
  }
  return(colors)
}



#' Generate colors of bar from https://doi.org/10.1016/j.cell.2020.03.048, Fig1A.
#' 
#' @param n number of colors. Max number is 44.
#' 
#' @return A color palette for plotting
#' 
#' @export
scPalette_anno_bar  <- function (n) 
{
  colorSpace <- c("#4493c2", "white", "#b11b2c")
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  }
  else {
    colors <- (grDevices::colorRampPalette(colorSpace))(n)
  }
  return(colors)
}
#' Generate colors of bar from https://doi.org/10.1016/j.cell.2020.03.048, Fig1A.
#' 
#' @param n number of colors. Max number is 44.
#' 
#' @return A color palette for plotting
#' 
#' @export
scPalette_anno_dark <- function (n) 
{
  colorSpace <- c("#caa37c", "#f1746c", "#e9a842", "#9fd6c7", 
                  "#ebaece", "#ba4a93", "#6bb2cd")
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  }
  else {
    colors <- (grDevices::colorRampPalette(colorSpace))(n)
  }
  return(colors)
}


#' Generate colors of bar from https://doi.org/10.1016/j.cell.2020.03.048, Fig1A.
#' 
#' @param n number of colors. Max number is 44.
#' 
#' @return A color palette for plotting
#' 
#' @export
scPalette_anno_light <- function (n) 
{
  colorSpace <- c("#ddc9b0", "#f6aba6", "#f2ca8d", "#9fd6c7", 
                  "#f4cee1", "#d592be", "#a3d1e1")
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  }
  else {
    colors <- (grDevices::colorRampPalette(colorSpace))(n)
  }
  return(colors)
}

#' Generate colors of heatmap bar.
#' 
#' @param color blue: red, white, blue; green: red, white, green.
#' 
#' @return A color palette for plotting
#' 
#' @export
scPalette_heatmap_bar <- function (color = c("blue", "green")[1]) 
{
  if (color == "blue") {
    colors <- (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, 
                                                                        name = "RdYlBu"))))(100)
    colors <- rev(colors)
  }
  if (color == "green") {
    colors <- c("#9E0242", "#BA2148", "#D6404E", "#E75947", 
                "#F57347", "#F99555", "#FDB567", "#FDCF7D", "#FFE695", 
                "#FDF7B0", "#F7FCB5", "#EBF7A0", "#D3ED9C", "#B4E1A2", 
                "#90D2A4", "#6EC5A4", "#51A9AF", "#358BBB", "#456CB0", 
                "#5E4FA2")
  }
  return(colors)
}


#' Generate colors inspired by R package: ggsci.
#' 
#' @param n Number of colors. Max number of colors per palette is 10.
#' 
#' @return A color palette for plotting
#' 
#' @examples 
#' scPalette_sci(10, palette = 'npg')
#' 
#' @export
scPalette_sci <- function (n, palette = c("npg", "aaas", "nejm", "lancet", "jama", 
                         "jco", "ucscgb", "d3", "locuszoom", "igv", "uchicago", "startrek", 
                         "tron", "futurama", "rickandmorty", "simpsons", "gsea", "material"), 
          ...) 
{
  choices = c("npg", "aaas", "nejm", "lancet", "jama", "jco", 
              "ucscgb", "d3", "locuszoom", "igv", "uchicago", "startrek", 
              "tron", "futurama", "rickandmorty", "simpsons", "gsea", 
              "material")
  all_pal = list(ggsci::pal_npg, ggsci::pal_aaas, ggsci::pal_nejm, 
                 ggsci::pal_lancet, ggsci::pal_jama, ggsci::pal_jco, ggsci::pal_ucscgb, 
                 ggsci::pal_d3, ggsci::pal_locuszoom, ggsci::pal_igv, 
                 ggsci::pal_uchicago, ggsci::pal_startrek, ggsci::pal_tron, 
                 ggsci::pal_futurama, ggsci::pal_rickandmorty, ggsci::pal_simpsons, 
                 ggsci::pal_gsea, ggsci::pal_material)
  names(all_pal) <- choices
  palette <- match.arg(arg = palette, choices = choices)
  if (!palette %in% choices) {
    stop("Invalid selection palette. Please choose one of palletes from ggsci.")
  }
  pal_fun <- all_pal[[palette]]
  maxnum <- suppressWarnings(sum(!is.na(pal_fun()(100))))
  if (n <= maxnum) {
    colors <- pal_fun(...)(n)
  }
  else {
    colors <- (grDevices::colorRampPalette(pal_fun(...)(maxnum)))(n)
  }
  return(colors)
}



#' Generate colors vector
#' 
#' @param n number of colors.
#' @param palette palette of colors. 3 palettes are provided in function: 1,2,3. Or input a vector by user.
#' @param color.ramps To extend color if n is larger than palette.
#' 
#' @export
color_generator <- function (n, palette = 1, color.ramps = F) {
  .Deprecated('angrycell::scPalette2')
  colorSpace1 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                   "#F29403", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4", 
                   "#222F75", "#1B9E77", "#B2DF8A", "#E3BE00", "#FB9A99", 
                   "#E7298A", "#910241", "#00CDD1", "#A6CEE3", "#CE1261", 
                   "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", "#B3DE69", 
                   "#8DD3C7", "#999999")
  colorSpace2 <- c("#00A087FF", "#4DBBD5FF", "#E64B35FF", "#3C5488FF", 
                   "#F38400", "#A1CAF1", "#BE0032", "#C2B280", "pink", "#008856", 
                   "#E68FAC", "#0067A5", "#604E97", "#F6A600", "#B3446C", 
                   "#DCD300", "#882D17", "#8DB600", "#654522", "#E25822", 
                   "#2B3D26", "#848482", "#008856", "#E68FAC", "#0067A5", 
                   "#604E97", "#F6A600", "#B3446C", "#DCD300", "#882D17")
  colorSpace3 <- c("#b589bc", "#ecb888", "#fc7440", "#a032cb", 
                   "#cd8ab1", "#e9d6eb", "#77a1ec", "#a3c4a5", "#eed785", 
                   "#d4b595", "#efe0e7", "#eebed6", "#f7d39b", "#cd7560", 
                   "#7fc28e", "#f9cabe", "#c76da8", "#fef6ae", "#6276b2", 
                   "#e69d8e", "#da7574", "#c5d9df", "#e16db7", "#908ebc", 
                   "#af88bb", "#dedbee", "#f07590", "#dfb9d5", "#b099b5", 
                   "#5394c3", "#b8a89f", "#917393", "#85c1d7", "#dbebf7", 
                   "#84d7f7", "#d4eef6", "#f5bcc9", "#5eb4f3", "#f7f4b7", 
                   "#6f89d3", "#50c2fb", "#ea94cd", "#d8b6f7", "#efc1ec")
  if (palette == 1 & n <= length(colorSpace1)) {
    colors <- colorSpace1[1:n]
  }
  else if (palette == 2 & n <= length(colorSpace2)) {
    colors <- colorSpace2[1:n]
  }
  else if (palette == 3 & n <= length(colorSpace3)) {
    colors <- colorSpace3[1:n]
  }
  else if (n <= length(colorSpace2)) {
    colors <- colorSpace2[1:n]
  }
  else if (n <= length(colorSpace3)) {
    colors <- colorSpace3[1:n]
  }
  else {
    colors <- (grDevices::colorRampPalette(colorSpace3))(n)
  }
  if (palette == 1 & color.ramps) {
    if (n <= length(colorSpace1)) {
      colors <- colorSpace1[1:n]
    }
    else {
      colors <- (grDevices::colorRampPalette(colorSpace1))(n)
    }
  }
  else if (palette == 2 & color.ramps) {
    if (n <= length(colorSpace2)) {
      colors <- colorSpace2[1:n]
    }
    else {
      colors <- (grDevices::colorRampPalette(colorSpace2))(n)
    }
  }
  else if (palette == 3 & color.ramps) {
    if (n <= length(colorSpace3)) {
      colors <- colorSpace3[1:n]
    }
    else {
      colors <- (grDevices::colorRampPalette(colorSpace3))(n)
    }
  }
  if (class(palette) != "numeric") {
    if (n <= length(palette)) {
      colors <- palette[1:n]
    }
    else {
      colors <- (grDevices::colorRampPalette(palette))(n)
    }
  }
  return(colors)
}


### other color palette
textblockcolor <- c('#ddc9b0','#f6aba6','#f2ca8d','#9fd6c7','#f4cee1','#d592be',
                    '#a3d1e1')

annocolor <- c('#caa37c','#f1746c','#e9a842','#9fd6c7','#ebaece','#ba4a93',
               '#6bb2cd')

barcolor <- c('#93c4dc','white','#f5b08f','#b93032','#b11b2c')
cellcolor <- c('#b589bc','#ecb888','#fc7440','#a032cb','#cd8ab1','#e9d6eb',
               '#77a1ec','#a3c4a5','#eed785','#d4b595','#efe0e7','#eebed6',
               '#f7d39b','#cd7560','#7fc28e','#f9cabe','#c76da8','#fef6ae',
               '#6276b2','#e69d8e','#da7574','#c5d9df','#e16db7','#908ebc',
               '#af88bb','#dedbee','#f07590','#dfb9d5','#b099b5','#5394c3',
               '#b8a89f','#917393','#85c1d7','#dbebf7','#84d7f7','#d4eef6',
               '#f5bcc9','#5eb4f3','#f7f4b7','#6f89d3','#50c2fb','#ea94cd',
               '#d8b6f7','#efc1ec')
heatmapcolor <- c('#4493c2','white','#b11b2c')
mycol <- c('#A50026','#D73027','#F46D43','#FDAE61','#FEE090','#FFFFBF',
           '#E0F3F8','#ABD9E9','#74ADD1','#4575B4','#313695')
