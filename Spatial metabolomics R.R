#!/opt/conda/bin/Rscript
#' @export
bubblescale <- function(data) {
  maxn <- max(data)
  minn <- min(data)
  if ((maxn - minn) == 0) {
    x <- maxn
    lab <- c(x)
    ran <- c(3, 3)
  } else if ((maxn - minn) == 1) {
    lab <- c(minn, minn + 1)
    ran <- c(3, 4)
  } else if ((maxn - minn) == 2) {
    lab <- c(minn, minn + 1, minn + 2)
    ran <- c(3, 5)
  } else if ((maxn - minn) %% 5 == 0 ||
             (maxn - minn + 1) %% 5 == 0) {
    x <- round((maxn - minn) / 5)
    lab <-
      c(minn, minn + x, minn + 2 * x, minn + 3 * x, minn + 4 * x, minn + 5 * x)
    ran <- c(1, 5)
  } else if ((maxn - minn) %% 4 == 0 ||
             (maxn - minn + 1) %% 4 == 0) {
    x <- round((maxn - minn) / 4)
    lab <-
      c(minn, minn + x, minn + 2 * x, minn + 3 * x, minn + 4 * x)
    ran <- c(1, 5)
  } else if ((maxn - minn) %% 3 == 0 ||
             (maxn - minn + 1) %% 3 == 0) {
    x <- round((maxn - minn) / 3)
    lab <- c(minn, minn + x, minn + 2 * x, minn + 3 * x)
    ran <- c(2, 5)
  } else {
    x <- floor((maxn - minn) / 4)
    lab <-
      c(minn, minn + x, minn + 2 * x, minn + 3 * x, minn + 4 * x)
    ran <- c(1, 5)
  }
  return(list(lab, ran))
}
#' 富集气泡图
#'
#' @param savepath 
#' @param type 
#' @param indata 
#' @param incompare 
#' @param filt 
#' @param grid 
#' @param number 
#' @param height 
#' @param width 
#' @param dpi 
#' @param fontfamily 
#' @param ...
#'
#' @export
bubble <- function(savepath = "./",
                   type = NULL,
                   indata = "None",
                   incompare = NULL,
                   filt = "T",
                   grid = "T",
                   number = 20,
                   imagetype = c("pdf", "png"),
                   height = 8,
                   width = 14,
                   dpi = 300,
                   fontfamily = "sans",
                   col = "brewer_spectra",
                   ...) {
  pacman::p_load(dplyr, ggplot2, stringr, Hmisc)
  mycol <- SelectColors(col)
  
  if (indata == "None") {
    comenrichpath <-
      paste0(savepath, getenrichtype(type)$filn, "/", incompare, "/")
    enrichfiles <-
      paste0(comenrichpath,
             dir(path = comenrichpath, pattern = "enrichment-*"))
    savepath = "./enrich/"
  } else{
    comenrichpath <- paste0(dirname(indata), "/")
    enrichfiles <- indata
  }
  typef <-
    ifelse(is.null(type), "enrichment", getenrichtype(type)$filn)
  fixed_col <- c("p-value", "ListHits", "Enrichment_score", "Term")
  for (enrichfile in enrichfiles) {
    tudd <-
      ifelse(grepl("-Total", enrichfile),
             "Total",
             ifelse(
               grepl("-Up", enrichfile),
               "Up",
               ifelse(grepl("-Down", enrichfile), "Down", "Total")
             ))
    endat <- readdata(enrichfile)
    file_col <- colnames(endat)
    stopname <- setdiff(fixed_col, file_col)
    if (length(stopname) > 0) {
      stop(paste0("NA", paste0(stopname, collapse = "、"), "Column or case names do not correspond"))
    }
    endat$`p-value`[endat$`p-value` == 0] <-
      min(endat$`p-value`[endat$`p-value` != 0])
    if (filt == "T") {
      if (grid == "T") {
        if (!"Category" %in% file_col) {
          stop("When selecting facets, the input table must have a Category column! Check if the Category column exists or if the names correspond")
        }
        endata_filt <-
          rbind(
            filter(endat, ListHits != 1, Category == "molecular_function")[1:5, ],
            filter(endat, ListHits != 1, Category == "cellular_component")[1:5, ],
            filter(endat, ListHits != 1, Category == "biological_process")[1:5, ]
          )
        endata <- endata_filt %>% arrange(Category, Enrichment_score)
      } else{
        endata_filt <- endat[order(endat$`p-value`), ][1:number, ]
        endata <- endata_filt %>% arrange(Enrichment_score)
      }
      na.omit(endata) -> endata
      savexls(endata,
              paste0(comenrichpath, typef, ".Bubble.", tudd, ".xls"))
    } else {
      endata <- endat
      na.omit(endata) -> endata
      endata <- endata[nrow(endata):1, ]
    }
    if (nrow(endata) > 2) {
      endata$Term <- capitalize(as.character(endata$Term))
      endata$Term <- factor(endata$Term, levels = endata$Term)
      if (grid == "T") {
        endata$Category <-
          str_replace_all(
            endata$Category,
            c(
              "biological_process" = "BP",
              "cellular_component" = "CC",
              "molecular_function" = "MF"
            )
          )
        
        p <- ggplot(endata, aes(Enrichment_score, Term)) +
          geom_point(aes(size = ListHits, color = endata$`p-value`)) +
          facet_grid(Category ~ ., scales = "free") +
          guides(size = guide_legend(order = 1)) +
          labs(
            x = "Enrichment_score",
            y = "",
            title = ifelse(
              is.null(incompare),
              paste0("Top ", typef, " Terms"),
              paste0(incompare, "(", tudd, ")  \n Top ", typef, " Terms")
            ),
            color = "P-value",
            size = "Count"
          ) +
          scale_colour_gradientn(colours = mycol) +
          theme_bw() +
          theme(aspect.ratio = 19 / 21) +
          theme(
            strip.background = element_rect(
              fill = "gray60",
              color = "gray60",
              linewidth = 1,
              linetype = "solid"
            )
          ) +
          theme(panel.border = element_rect(
            color = "black",
            linewidth = 1,
            fill = NA
          )) +
          scale_size_continuous(breaks = bubblescale(endata$ListHits)[[1]],
                                range = bubblescale(endata$ListHits)[[2]]) +
          scale_y_discrete(
            labels = function(x)
              ifelse(nchar(x) > 60, paste0(substr(x, 1, 60), "..."), x)
          ) +
          theme(plot.title = element_text(size = 18, hjust = 0.5)) +
          theme(plot.margin = unit(c(2, 0, 2, 0), "lines")) +
          theme(
            axis.text.x = element_text(size = 15, color = "black", margin = margin(t = 15, b = 0)),
            strip.text = element_text(size = 15, color = "black"),
            axis.text.y = element_text(size = 15, color = "black"),
            axis.title = element_text(size = 18, color = "black"),
            axis.title.x = element_text(margin = margin(t = 20, b = 0)),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14)
          )
        
      } else{
        p <- ggplot(endata, aes(RichFactor, Term)) +
          geom_point(aes(size = ListHits, color = endata$`p-value`)) +
          scale_y_discrete(
            labels = function(x)
              ifelse(nchar(x) > 60, paste0(substr(x, 1, 60), "..."), x),
            position = "right"
          ) +
          theme_bw() +
          theme(panel.border = element_rect(
            color = "black",
            linewidth = 1,
            fill = NA
          )) +
          theme(axis.line = element_blank(),
                strip.background = element_blank()) +
          theme(aspect.ratio = 16 / 9) +
          theme(
            axis.text.x = element_text(size = 12, color = "black", margin = margin(t = 10, b = 0)),
            axis.text.y = element_text(size = 15, color = "black"),
            axis.title = element_text(size = 16, color = "black"),
            axis.title.x = element_text(margin = margin(t = 15, b = 0)),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            legend.position = c(1.5, -0.05),
            legend.direction = "horizontal",
            legend.box = "horizontal",
            legend.box.just = "right",
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
            legend.background = element_rect(fill = "white", color = NA),
            legend.spacing.x = unit(1, "cm"),
            plot.margin = unit(c(1, 8, 1, 1), "lines")
          ) +
          scale_colour_gradientn(colours = mycol) +
          scale_size_continuous(breaks = bubblescale(endata$ListHits)[[1]],
                                range = bubblescale(endata$ListHits)[[2]]) +
          guides(
            size = guide_legend(
              order = 1, 
              direction = "horizontal", 
              title.position = "bottom",
              title.hjust = 0.5,
              label.position = "top",
              keywidth = unit(1, "lines"),
              keyheight = unit(1, "lines")
            ),
            color = guide_colorbar(
              direction = "horizontal", 
              title.position = "bottom",
              title.hjust = 0.5,
              label.position = "top",
              barwidth = unit(5, "cm"),
              barheight = unit(0.8, "lines")
            )
          ) +
          labs(
            x = "RichFactor",
            y = '',
            title = NULL,
            color = "P-value",
            size = "Count"
          )
        
      }
      ggplotsave(
        plot = p,
        savepath = comenrichpath,
        mapname = paste0(typef, ".Bubble.", tudd),
        width = width,
        height = height,
        imagetype = imagetype,
        family = fontfamily,
        ...
      )
      
      
    } else{
      warning("The number of terms is less than 3 and no enrichment bubble chart is provided")
      savetxt(
        data = "The number of terms is less than 3 and no enrichment bubble chart is provided",
        filename = paste0(comenrichpath, "/说明.txt"),
        append = T
      )
    }
  }
}

if (this.path::is.main() & !is.null(whereami::thisfile_r())) {
  options(warn = -1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # parameters
  parser$add_argument("-i",
                      "--imagetype",
                      default = c("pdf", "png"),
                      help = "Image format")
  parser$add_argument("-fa", "--fontfamily", default = "sans", help = "Font, default is Arial")
  parser$add_argument("-wi",
                      "--width",
                      default = 14,
                      type = "double",
                      help = "Image width")
  parser$add_argument(
    "-he",
    "--height",
    default = 8,
    type = "double",
    help = "Image height"
  )
  parser$add_argument("-d",
                      "--dpi",
                      default = 300,
                      type = "double",
                      help = "resolution")
  parser$add_argument("-col", "--col", default = "RdYlGn" , help = "Color scheme of bubble diagram")
  
  parser$add_argument(
    "-t",
    "--type",
    type = "character",
    default = NULL,
    help = "enrich database class,such as G/K/W/R/I",
    metavar = "character"
  )
  parser$add_argument("-ic", "--incompare", default = NULL, help = "Compare group names, default NULL")
  parser$add_argument("-id", "--indata", default = "None", help = "Drawing data table file -LRB-include path) , default None is not entered")
  parser$add_argument("-f",
                      "--filt",
                      default = "T",
                      type = "character",
                      help = "Whether the data needs to be filtered, default T")
  parser$add_argument("-g",
                      "--grid",
                      type = "character",
                      default = "F",
                      help = "If facets are required, defaults to F")
  parser$add_argument("-n", "--number", default = 20, help = "top number each category,default 20")
  parser$add_argument("-s", "--savepath", default = "./", help = "Enriched result save path, default. /, the default in the process. /enrich/")
  
  args <- parser$parse_args()
  
  bubbleplot <- do.call(bubble, args = args)
}

install.packages("BiocManager")
BiocManager::install("Cardinal")
library(Cardinal)

# Loading MSImagingArrays data
msa_file <- CardinalIO::exampleImzMLFile("processed")
msa <- readMSIData(msa_file)
msa

# Loading MSImagingExperiments data
mse_files <- CardinalIO::exampleImzMLFile("continuous")
mse <- readMSIData(mse_files)
mse

# Loading in house data
final <- readMSIData("final/final-neg.imzML")
final
pre <- readMSIData("pre/pre-neg.imzML")
pre

#Building MSImagingExperiment
set.seed(2020, kind="L'Ecuyer-CMRG")
mse <- simulateImage(preset=6, dim=c(40,40), baseline=0.5)
mse


#Convert to MSImagingArrays
msa <- convertMSImagingExperiment2Arrays(mse)

imageData(final)

# Accessing feature data
featureData(pre)
fData(final)

# Accessing spectra data
spectraData(pre)
spectra(final)

# Accessing pixel data
pData(pre)
pixelData(final)

plot(final)
plot(final, i=c(236,956))

plot(final, coord=list(x=12, y=7))

plot(final, xlim=c(200,400))

image(final)
image(final,i=2)
image(pre,mz = 143.1088)

image(final, mz = 143.1088, tolerance=1, units = "mz")

# Smooth
image(final,mz = 143.1088, smooth = "gaussian")
image(final,mz = 143.1088, smooth = "guide")

# Enhance
image(final, mz = 143.1088, enhance = "adaptive")
image(final, mz = 143.1088, enhance = "histogram")

cap <- selectROI(final,mz = 143.1088)

# Subset
final[1:50]
final[1:50,1:50]

subset(final, 70 < mz & mz < 100, x <= 10 & y <= 10)

# Make factors
mid <- readRDS("DESCK6-Mid1-neg.rds")
cap <- readRDS("DESCK6-Cap1-neg.rds")
other <- !mid & !cap

final$area <- makeFactor(mid=mid,cap=cap,other=other)
pData(final)

final_cap <- final[,final$area == "cap"]
final_mid <- final[,final$area == "mid"]

# calculate mean spectrum
final <- summarizeFeatures(final, stat = "mean")
fData(final)
plot(final,"mean")

# calculate TIC
final <- summarizePixels(final, stat=c(TIC="sum"))
pData(final)
image(final,"TIC")

final <- summarizeFeatures(final,stat = "mean", groups = final$area)
fData(final)
plot(final,c("mid.mean","cap.mean"), superpose = T)

# Normalization
msa_pre <- normalize(msa, method="tic")

# Smoothing
msa_smooth <- smooth(msa_pre, method="gaussian")
plot(msa_smooth,coord=list(x=9, y=8),
     xlim=c(1150,1350), linewidth=2)

# Peak picking
p1 <- peakPick(msa_pre, method="diff", SNR=3)
plot(p1,coord=list(x=16, y=16), linewidth=2)

# Peak alignment
msa_peaks <- peakPick(msa_pre, method="filter", SNR=3)
msa_peaks
mse_peaks <- peakAlign(msa_peaks, tolerance=200, units="ppm")
mse_peaks

# Peak filtering
mse_filt <- subsetFeatures(mse_peaks, freq > 0.1)
mse_filt

mse_filt_sum <- summarizeFeatures(mse_filt)

plot(mse_filt_sum, "mean", xlab="m/z", ylab="Intensity",
     linewidth=2, annPeaks=10)

#!/opt/conda/bin/Rscript
#' Multiple sets of volcanic maps are calculated FC
#' @export
Mulgroups_FC <- function(filename="data matrix.xlsx"){
  
  library(readxl)
  library(openxlsx)
  library(tidyr)
  
  data_mat <- read_excel(filename, sheet = "data matrix")
  group <- read_excel(filename, sheet = "group")
  data <- data_mat[, colnames(data_mat) %in% group$Sample]
  groups <- unique(group$Group)
  
  # 3. Calculate the FC for each group
  result <- NULL
  for (grp in groups) {
    # Gets the sample data corresponding to the current group
    group_samples <- group$Sample[group$Group == grp]
    
    group_data <- data[, group_samples]
    other_group_data <- data[, !colnames(data) %in% group_samples]
    # Calculates the fold change of the merged data of the current group and other groups
    mean_data1 <- apply(group_data, 1, mean, na.rm = TRUE)
    mean_data2 <- apply(other_group_data, 1, mean, na.rm = TRUE)
    fc <- log2(mean_data1 / mean_data2)
    # Construct the result data box
    df <- data.frame(
      mz = data_mat$mz,
      Group = grp,
      log2FC = fc
    )
    
    # Consolidated results
    if (is.null(result)) {
      result <- df
    } else {
      result <- rbind(result, df)
    }
  }
  
  write.xlsx(result, "volcano_data.xlsx", rowNames = F)
}


#' Multiple sets of volcanic maps
#' 
#' MulgroupsVolcano
#' 
#' @param filename 
#' @param savepath 
#' @param savename 
#' @param imagetype 
#' @param width 
#' @param height 
#' @param family 
#' @param myMarkers 
#' @param order.by 
#' @param log2FC.cutoff 
#' @param topGeneN 
#' @param back.col 
#' @param pSize 
#' @param aesCol 
#' @param legend.position 
#' @param base_size 
#' @param tile.col 
#' @param classfile 
#' @param group.order 
#' @param polar
#' @param expand 
#' @param flip
#' @param tilemarker 
#' @param grouptextsize 
#' @param titlesize 
#' @param linesize 
#' @param ticksize
#' @param legendsize 
#' @export
MulgroupsVolcano<-function(filename="volcano_data.xlsx",
                           savepath = "./", 
                           savename="Volcano",    
                           imagetype=c('png', 'pdf'),
                           width=10,
                           height=10,
                           family="sans",
                           myMarkers = NULL,
                           order.by = "log2FC", 
                           log2FC.cutoff = 0.25,
                           topGeneN = 5,   
                           back.col = 'grey95',
                           pSize = 1,
                           aesCol = c('#0099CC','#CC3333'),
                           legend.position = c(0.9,0.9),
                           base_size = 14,
                           classfile = "classtype.xlsx",
                           tile.col = stylefun_group(classfile = classfile,styletype = "fill",fillpalette = "procol"), 
                           group.order = NULL,
                           polar = FALSE,
                           expand = c(-1,1),
                           flip = FALSE,
                           tilemarker = NULL,
                           grouptextsize=6,
                           textsize=10,
                           titlesize=13,
                           linesize=1,
                           ticksize=1,
                           legendsize=13,
                           ...){
  library(rlang)
  library("magrittr")
  library("ggplot2")
  
  Mulgroups_FC(filename=filename)  
  diff.marker <- readdata(filename="volcano_data.xlsx")
  
  
  
  # assign type
  diff.marker <- diff.marker %>%
    dplyr::mutate(type = ifelse(log2FC >= log2FC.cutoff,"sigUp","sigDown"))
  
  # cluster orders
  if(!is.null(group.order)){diff.marker$Group <- factor(diff.marker$Group, levels = group.order)}
  
  # get background cols
  purrr::map_df(unique(diff.marker$Group),function(x){
    tmp <- diff.marker %>% dplyr::filter(Group == x)
    new.tmp <- data.frame(Group = x,
                          min = min(tmp$log2FC) - 0.2,
                          max = max(tmp$log2FC) + 0.2)
    return(new.tmp)
  }) -> back.data
  
  # get top gene
  top.marker.tmp <- diff.marker %>%  dplyr::group_by(Group)
  
  # order
  top.marker.max <- top.marker.tmp %>%
    dplyr::slice_max(n = topGeneN,order_by = get(order.by))
  
  top.marker.min <- top.marker.tmp %>%
    dplyr::slice_min(n = topGeneN,order_by = get(order.by))
  
  # combine
  top.marker <- rbind(top.marker.max,top.marker.min)
  
  
  # tilemarker gene
  tilemarkerdata <- top.marker %>% dplyr::filter(mz %in% tilemarker)
  
  # whether supply own genes
  if(!is.null(myMarkers)){
    top.marker <- diff.marker %>% dplyr::filter(mz %in% myMarkers)
  }else{
    top.marker <- top.marker
  }
  
  # plot
  p1 <- ggplot2::ggplot(diff.marker, ggplot2::aes(x = Group,y = log2FC)) +
    # add back cols
    ggplot2::geom_col(data = back.data,
                      ggplot2::aes(x = Group,y = min),fill = back.col) +
    ggplot2::geom_col(data = back.data,
                      ggplot2::aes(x = Group,y = max),fill = back.col)+
    theme(
      
      axis.text=element_text(size=textsize,color="black",face="bold"),
      axis.title=element_text(size=titlesize,color="black",face="bold"),
      axis.line.y= element_line(size = linesize),
      axis.ticks.x= element_line(size = ticksize)
    )
  
  # color type
  p2 <- p1 +
    # add point
    ggplot2::geom_jitter(ggplot2::aes(color = type),size = pSize) +
    ggplot2::scale_color_manual(values = c("sigDown" = aesCol[1],"sigUp" = aesCol[2]))
  
  
  # theme details
  p3 <- p2 +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = legend.position,
                   legend.title = ggplot2::element_blank(),
                   legend.background = ggplot2::element_blank()) +
    ggplot2::xlab("Group") + ggplot2::ylab('Log2(FoldChange)') +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5), ncol = 1, reverse = T))
  
  # if (min(back.data$min) > 0){
  #   p3.1 <- p3+ggbreak::scale_y_break(c(10,(max(back.data$min)-10)),
  #                          space = 0.1)
  #   
  #   
  #   gg.gap(plot =  p3,
  #          segments = c(0,(max(back.data$min)-10)),
  #          tick_width = 10,
  #          ylim = c(0,max(back.data$max)))
  #   
  # }
  
  
  
  
  # add tile
  p4 <- p3 +
    ggplot2::geom_tile(ggplot2::aes(x = Group,y = 0,fill = Group),
                       color = 'white',
                       height = log2FC.cutoff*2,
                       alpha = 0.3,
                       show.legend = F) +
    ggplot2::scale_fill_manual(values = tile.col) +
    # add gene label
    ggrepel::geom_text_repel(data = tilemarkerdata,
                             ggplot2::aes(x = Group,y = log2FC,label = mz, fontface = "bold"),
                             max.overlaps = 50)
  
  # whether coord_plolar
  if(polar == TRUE){
    p5 <- p4 +
      geomtextpath::geom_textpath(ggplot2::aes(x = Group,y = 0,label = Group)) +
      ggplot2::scale_y_continuous(n.breaks = 6,
                                  expand = ggplot2::expansion(mult = expand)) +
      ggplot2::theme_void(base_size = base_size) +
      ggplot2::theme(legend.position = legend.position,
                     legend.title = ggplot2::element_blank()) +
      ggplot2::coord_polar(clip = 'off',theta = 'x')
  }else{
    # whether flip plot
    if(flip == TRUE){
      p5 <- p4 +
        ggplot2::scale_y_continuous(n.breaks = 6) +
        ggplot2::geom_label(ggplot2::aes(x = Group,y = 0,label = Group), color="white") +
        ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank()) +
        ggplot2::coord_flip()
    }else{
      p5 <- p4 +
        ggplot2::scale_y_continuous(n.breaks = 6) +
        ggplot2::geom_text(ggplot2::aes(x = Group,y = 0,label = Group), color="white",size=grouptextsize) +
        ggplot2::theme(axis.line.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.text=element_text(size=textsize,color="black",face="bold"),
                       axis.title=element_text(size=titlesize,color="black",face="bold"),
                       axis.line.y= element_line(size = linesize),
                       axis.ticks.y= element_line(size = ticksize),
                       legend.text=element_text(size=legendsize,color="black",face="bold"))
    }
  }
  
  
  
  ggplotsave(plot = p5,limitsize = FALSE,
             savepath = savepath,
             mapname = savename,
             imagetype = imagetype,
             family = family,
             width = width,
             height = height)
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  suppressMessages(library("rlang"))
  suppressMessages(library("magrittr"))
  suppressMessages(library("ggplot2"))
  
  parser <- ArgumentParser()
  parser$add_argument("-f","--filename",default="data matrix.xlsx",help = "data matrix.xlsx:Contains a data matrix and a grouped sheet")
  parser$add_argument("-sp","--savepath",default = "./",help = "Save Results")
  parser$add_argument("-sm","--savename", default = "Volcano", help = "Name of drawing")
  parser$add_argument("-it","--imagetype",default = c("jpg","pdf"), help = "Image format",nargs="+",choices = c("jpg","tiff","png","pdf"))
  parser$add_argument("-fa","--family",default = "sans", help = "font")
  parser$add_argument("-wi","--width",default = 10, type= "double",help = "Image width")
  parser$add_argument("-he","--height",default = 10, type= "double",help = "Image height")
  parser$add_argument("-cf","--classfile",default = "classtype.xlsx", help = "Group color templates")
  parser$add_argument("-mM","--myMarkers",default=NULL,nargs="+",help="Select only specific data")
  parser$add_argument("-ob","--order.by ",default="log2FC",help="Data sorting method")
  parser$add_argument("-lc","--log2FC.cutoff",default=0.25,help="Criterion of significant difference")
  parser$add_argument("-topn","--topGeneN",default=5,help="top")
  parser$add_argument("-bc","--back.col",default="grey95",help="Background color of volcano image")
  parser$add_argument("-ps","--pSize",default=1,help="The size of the dot")
  parser$add_argument("-ac","--aesCol",default=c('#0099CC','#CC3333'),nargs="+",help="Up and down color, first down, then up")
  parser$add_argument("-lp","--legend.position",default=c(0.9,0.9),nargs="+",help="Legend location")
  parser$add_argument("-bs","--base_size",default=14,help="Base font size for drawing")
  parser$add_argument("-go","--group.order",default=NULL,nargs="+",help="Grouping and sorting of volcano map")
  parser$add_argument("-ep","--expand",default=c(-1,1),nargs="+",help="Coefficient of expansion of coordinate axis")
  parser$add_argument("-fl","--flip",default=FALSE,help="Whether to rotate the coordinate system")
  parser$add_argument("-po","--polar",default=FALSE,help="Whether or not to convert to polar coordinates")
  parser$add_argument("-tm","--tilemarker",default=NULL,nargs="+",help="The data to be annotated")
  parser$add_argument("-gts","--grouptextsize",default=6,help="Group font sizes")
  parser$add_argument("-tes","--textsize",default=10,help="Axis label font size")
  parser$add_argument("-tits","--titlesize",default=13,help="Axis heading font size")
  parser$add_argument("-lis","--linesize",default=1,help="Thickness of the coordinate axis")
  parser$add_argument("-tics","--ticksize",default=1,help="Thickness of the graduation line")
  parser$add_argument("-les","--legendsize",default=13,help="Legend font size")
  
  
  args <- parse_args(parser)
  
  
  result <- do.call(what = MulgroupsVolcano, args = args)
  
}
