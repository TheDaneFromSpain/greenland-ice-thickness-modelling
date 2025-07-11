my.create_uniform_layout <- function(n, orientation = "landscape", nrow = NULL, ncol = NULL, numbering = seq(1, dim1*dim2), ...){
  if (n>1){
    if (!is.null(nrow)){
      dim1 <- nrow
    }else{
      dim1 <- which.min((1:n - sqrt(n))^2)
    }
    if (!is.null(ncol)){
      dim2 <- ncol
      dim1 <- ceiling(n/dim2)
    }else{
      dim2 <- ceiling(n/dim1)
    }
    if (orientation == "portrait" | dim1 >= dim2){
      layout <- layout(matrix(numbering, max(dim1, dim2), min(dim1, dim2), ...))
    }else if (orientation == "landscape" | dim1 < dim2){
      layout <- layout(matrix(numbering, min(dim1, dim2), max(dim1, dim2), ...))
    }else{
      message("Orientation must be landscape or portrait")
    }
    layout.show(layout)
  }
}

my.get_consistent_colors <- function(range_matrix,  # rows: min, max; cols: data1, data2, ...
                                     amount = 256,  # number of colors
                                     palette = "viridis",  # hcl palette to use
                                     reverse_palette = FALSE
){ 
  n <- ncol(range_matrix)
  m <- min(range_matrix); M <- max(range_matrix)
  len <- M - m
  full_palette = hcl.colors(amount, palette, rev = reverse_palette)
  
  colors = list()
  for (i in 1:n){ 
    colors[[i]] <- full_palette[round(amount * (range_matrix[1, i] - m) / len):round(amount * (1 - (M - range_matrix[2, i]) / len))]
  }
  
  names(colors) <- colnames(range_matrix)
  return(colors)
}

my.plot_histogram_boxplots <- function(data.df_data.numlist, # data frame with data in columns OR list of numeric vectors
                                       global_x = NULL, # xlab if all x use the same variable
                                       breaks = "Sturges", # "Scott"
                                       filename = "histogram_boxplot", # output png file name
                                       resolution = c(1920, 1080), # image resolution: width, height
                                       ...
){
  png(paste0(filename, ".png"), width = resolution[1], height = resolution[2])
  my.create_uniform_layout(length(names(data.df_data.numlist)), ...)
  
  for (name in names(data.df_data.numlist)){
    data <- data.df_data.numlist[[name]]
    
    if (is.null(global_x)){
      hist(data, breaks = breaks, main = "", ylab = "Frequency", xlab = gsub("_", " ", name), cex.lab = 1.5)
    }else{
      hist(data, breaks = breaks, main = gsub("_", " ", name), ylab = "Frequency", xlab = global_x, cex.lab = 1.5)
    }
    par(new = TRUE)
    boxplot(data, main = "", horizontal = T, axes = F, col = rgb(0, 0.8, 1, 0.5), cex = 2)
  }
  
  dev.off()
}

my.plot_scatterplots <- function(variable.num, # vector of observed variable of interest values
                                 variable_name = "Variable of interest",
                                 features.df, # dataframe with columns being feature values at observed variable locations
                                 filename = paste0("scatterplots_", variable_name), # output png file name
                                 resolution = c(1920, 1080), # image resolution: width, height
                                 palette = "rocket", # hcl palette to use
                                 amount = 256, # number of colours to use
                                 reverse_palette = FALSE,
                                 pointsize = 2,
                                 spread = 50, # controls the spread of colour in the scatterplot
                                 ...
){
  variable_range <- diff(range(variable.num))
  
  png(paste0(filename, ".png"), width = resolution[1], height = resolution[2], pointsize = 12*pointsize)
  my.create_uniform_layout(ncol(features.df), ...)
  
  for (name in colnames(features.df)){
    feature.num <- features.df[, name]
    
    feature_range <- diff(range(feature.num))
    correlation <- round(cor(feature.num, variable.num), digits = 3)
    
    scattercolours <- densCols(feature.num, variable.num, nbin = amount, 
                               colramp = colorRampPalette(hcl.colors(amount, palette, rev = reverse_palette)), 
                               bandwidth = c(feature_range / spread, variable_range / spread))
    
    plot(feature.num, variable.num, main = paste0("Correlación: ", format(correlation, nsmall = 3)),
         xlab = name, ylab = variable_name, cex.lab = 1.5, col = scattercolours)
  }
  
  dev.off()
}

my.plot_crossplot_matrix <- function(features.df, # dataframe with columns being feature values at observed variable locations
                                     filename = "feature_crossplot_matrix", # output png file name
                                     resolution = c(1920, 1080), # image resolution: width, height
                                     palette = "rocket", # hcl palette to use
                                     amount = 256, # number of colours to use
                                     reverse_palette = FALSE,
                                     pointsize = 2,
                                     spread = 50 # controls the spread of colour in the scatterplot
){
  png(paste0(filename, ".png"), width = resolution[1], height = resolution[2], pointsize = 12*pointsize)

  correlation_panel <- function(x, y){
    # usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y), digits = 3)
    txt <- paste0("corr = ", format(r, nsmall = 3))
    cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * sqrt(abs(r) + .5))
  }
  
  scatterplot_panel <- function(x, y){
    
    scattercolours <- densCols(x, y, nbin = amount, 
                               colramp = colorRampPalette(hcl.colors(amount, palette, rev = reverse_palette)), 
                               bandwidth = c(diff(range(x)) / spread, diff(range(y)) / spread))
    
    points(x, y, cex.lab = 1.5, col = scattercolours)
  }
  
  pairs(features.df, 
        lower.panel = correlation_panel,
        upper.panel = scatterplot_panel)
  
  dev.off()
}

my.plot_ascending_boxplot_variance <- function(data.num,
                                               data_name = "data",
                                               nsplits = 20,
                                               filename = paste0("boxplot_variance_", data_name), # output png file name
                                               resolution = c(1920, 1080) # image resolution: width, height
){
  png(paste0(filename, ".png"), width = resolution[1], height = resolution[2])
  
  split_data <- split(sort(data.num), ceiling(seq_along(data.num)/ceiling(length(data.num)/nsplits)))
  boxplot(split_data, axes = F, horizontal = TRUE,
          main = paste0(data_name, " boxplots with variance"), col = rgb(0, 0.8, 1, 0.5))
  par(new=TRUE)
  # plot(vars, type="b", pch = 19, cex = .5, ann = F, axes = F, xgap.axis = )
  plot(diff(range(data.num)) * sapply(split_data, median) / max(unlist(split_data)), sapply(split_data, var),
       type = "b", col = "red", pch = 19, yaxt = "n", ylab = "", xlab = data_name, xlim = c(0, max(unlist(split_data))))
  axis(2, col.ticks = "red", col.axis = "red")
  mtext("Variance", side = 2, line = 3, col = "red")
  par(new=FALSE)
  
  dev.off()
}

my.plot_terra_maps <- function(spat.list, # list of spatvectors or rasters to plot
                               study_variable = "raster", # name of variable to be plotted
                               background.rast = NULL, # background raster to display
                               max_points = Inf, # max number of points of spatvect to plot per map
                               palette = "Temps",
                               amount = 256,
                               consistent_colors = TRUE,
                               reverse_palette = FALSE,
                               cex = .1,
                               pointsize = 2,
                               filename = paste0("all_", study_variable, "_maps"), # output png file name
                               resolution = c(1920, 1080), # image resolution: width, height
                               ...
){
  png(paste0(filename, ".png"), width = resolution[1], height = resolution[2], pointsize = 12*pointsize)
  my.create_uniform_layout(length(spat.list), ...)
  
  if (length(dim(spat.list[[1]])) == 3){
    
    sapply(names(spat.list), function(name){ plot(spat.list[[name]], main = gsub("_", " ", name), 
                                                  legend = FALSE, # for personalized legend
                                                 col = hcl.colors(amount, palette = palette, rev = reverse_palette));
           legend("bottomright", legend = c("Ocean", "Ice-free land     ", "Grounded ice     ", # for personalized legend
                                            "Floating ice", "Canada"), fill = hcl.colors(5, palette = palette, rev = reverse_palette)) # for personalized legend
      # legend("bottomright", legend = c("Other", "Altitude", "Mass conserv.", "Kriging & Fluid     ") # for personalized legend
      #                                  , fill = hcl.colors(4, palette = palette, rev = reverse_palette)) # for personalized legend
      })
    
  }else{
    
    range_matrix <- matrix(unlist(lapply(spat.list, function(vect) range(vect[[study_variable]][,1]))), nrow = 2)
    if (!is.null(palette)){
      if (consistent_colors){
        col_palette <- my.get_consistent_colors(range_matrix = range_matrix, amount = amount, palette = palette, reverse_palette = reverse_palette)
        names(col_palette) <- names(spat.list)
      }else{
        col_palette <- lapply(1:length(spat.list), function(x) hcl.colors(amount, palette = palette, rev = reverse_palette))
        names(col_palette) <- names(spat.list)
      }
    }else{
      col_palette <- NULL
    }
    
    if (is.null(background.rast)){
      sapply(names(spat.list), function(name){
        plot(sample(spat.list[[name]], min(length(spat.list[[name]]), max_points)),
             study_variable, main = gsub("_", " ", name), type = "continuous", cex = cex, col = col_palette[[name]]) })
    }else{
      sapply(names(spat.list), function(name){
        plot(background.rast, col = gray(seq(.1, .9, by = .01)), main = gsub("_", " ", name), legend = FALSE)
        plot(sample(spat.list[[name]], min(length(spat.list[[name]]), max_points)),
             study_variable, type = "continuous", cex = cex, col = col_palette[[name]], add = TRUE) })
    }
  }
  
  dev.off()
}

my.plot_terra_maps_histograms_boxplots <- function(spat.list, # list of spatvectors OR rasters to plot
                                                   study_variable = "raster", # name of variable to be plotted
                                                   background.rast = NULL, # background raster to display
                                                   max_points = Inf, # max number of points of spatvect to plot per map
                                                   palette = "Temps",
                                                   amount = 256,
                                                   consistent_colors = TRUE,
                                                   reverse_palette = FALSE,
                                                   global_x = NULL, # xlab if all x use the same variable
                                                   breaks = "Sturges", # "Scott"
                                                   cex = .1,
                                                   pointsize = 2,
                                                   filename = paste0("all_", study_variable, "_maps_histograms_boxplots"), # output png file name
                                                   resolution = c(1920, 1080), # image resolution: width, height
                                                   ...
){
  png(paste0(filename, ".png"), width = resolution[1], height = resolution[2], pointsize = 12*pointsize)
  my.create_uniform_layout(length(spat.list)*2, ...)
  
  if (length(dim(spat.list[[1]])) == 3){
    
    for (name in names(spat.list)){
      
      plot(spat.list[[name]], main = gsub("_", " ", name), col = hcl.colors(amount, palette = palette, rev = reverse_palette))
      
      data <- as.data.frame(spat.list[[name]])[,1]

      hist(data, breaks = breaks, main = "", ylab = "Frequency", xlab = gsub("_", " ", name), cex.lab = 1.5)
      par(new = TRUE)
      boxplot(data, main = "", horizontal = T, axes = F, col = rgb(0, 0.8, 1, 0.5), cex = 2)
    }
    
  }else{
    
    range_matrix <- matrix(unlist(lapply(spat.list, function(vect) range(vect[[study_variable]][,1]))), nrow = 2)
    if (!is.null(palette)){
      if (consistent_colors){
        col_palette <- my.get_consistent_colors(range_matrix = range_matrix, amount = amount, palette = palette, reverse_palette = reverse_palette)
        names(col_palette) <- names(spat.list)
      }else{
        col_palette <- lapply(1:length(spat.list), function(x) hcl.colors(amount, palette = palette, rev = reverse_palette))
        names(col_palette) <- names(spat.list)
      }
    }else{
      col_palette <- NULL
    }
    
    if (is.null(background.rast)){
      for (name in names(spat.list)){
        
        vect <- sample(spat.list[[name]], min(length(spat.list[[name]]), max_points))
        
        plot(vect, study_variable, main = gsub("_", " ", name), type = "continuous", cex = cex, col = col_palette[[name]])
      
        data <- as.data.frame(vect)[,study_variable]
        
        if (is.null(global_x)){
          hist(data, breaks = breaks, main = "", ylab = "Frequency", xlab = gsub("_", " ", name), cex.lab = 1.5)
        }else{
          hist(data, breaks = breaks, main = gsub("_", " ", name), ylab = "Frequency", xlab = global_x, cex.lab = 1.5)
        }
        par(new = TRUE)
        boxplot(data, main = "", horizontal = T, axes = F, col = rgb(0, 0.8, 1, 0.5), cex = 2)
      }
    }else{
      for (name in names(spat.list)){
      
        vect <- sample(spat.list[[name]], min(length(spat.list[[name]]), max_points))
        
        plot(background.rast, col = gray(seq(.1, .9, by = .01)), main = gsub("_", " ", name), legend = FALSE)
        plot(vect, study_variable, type = "continuous", cex = cex, col = col_palette[[name]], add = TRUE)
      
        data <- as.data.frame(vect)[,study_variable]
      
        if (is.null(global_x)){
          hist(data, breaks = breaks, main = "", ylab = "Frequency", xlab = gsub("_", " ", name), cex.lab = 1.5)
        }else{
          hist(data, breaks = breaks, main = gsub("_", " ", name), ylab = "Frequency", xlab = global_x, cex.lab = 1.5)
        }
        par(new = TRUE)
        boxplot(data, main = "", horizontal = T, axes = F, col = rgb(0, 0.8, 1, 0.5), cex = 2)
      }
    }
  }
  dev.off()
}