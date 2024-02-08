

#' this function is called by others to compute the localizations of grid points
#' @param spatial_coords a dataframe with cell/spot as row and axis (X, Y) as column
#' @param n.grids.scale number of grids to use equals to (total number of cells)/n.grids.scales
#' @return returns a dataframe with grid points as rows and coordinates as columns


get.grid.coords <- function(spatial_coords, n.grids.scale = 5){
  n.grids <- floor(nrow(spatial_coords)/n.grids.scale)
  cl <- balanced_clustering(spatial_coords, K = n.grids)
  centroid.coords <- aggregate(spatial_coords, list(cl), mean)[, -1]
  return(centroid.coords)
}

#' from an excellent post: https://www.r-bloggers.com/2013/05/pairwise-distances-in-r/
#' this function is called by other functions to quickly compute the distance between
#' cells to grid points, or between grid points
#'
#' @param A matrix
#' @param B matrix
#' @return returns pairwise-distances
#'
spa_vectorized_pdist <- function(A,B){
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))

  m = nrow(A)
  n = nrow(B)

  tmp = matrix(rep(an, n), nrow=m)
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  return(sqrt( tmp - 2 * tcrossprod(A,B) ))
}

#' calculates local gradients at each grid point
#'
#' @param grids output of get.grid.coords
#' @param dist.to.grid distance between each cell/spot to each grid point
#' @param latent.embeddings embeddings used to calculate gradient with
#' @param spatial_coords a dataframe with cell/spot as row and axis (X, Y) as column
#' @param n.cells.per.meta the size of local group; it is suggested to use 20-50 as a start
#' @return returns local linear gradient information
#'
calc.local.linear <- function(grids, dist.to.grid, latent.embeddings, spatial_coords, n.cells.per.meta){
  cell.list <- list()
  assignment.list <- list()
  coeff.list <- list()
  # coeff.list.latent <- list()
  rsquared.list <- list()
  full.res.list <- list()
  colnames(spatial_coords) <- c("X", "Y")
  for (i in 1:ncol(dist.to.grid)){
    o <- order(dist.to.grid[, i])[1:n.cells.per.meta]
    cells <- rownames(dist.to.grid)[o]
    cell.list[[i]] <- cells
    local.spatial_coords <- spatial_coords[cells, ]
    local.spatial_coords <- scale(local.spatial_coords)
    local.latent.embeddings <- latent.embeddings[cells, ]
    calc.df <- as.data.frame(cbind(local.latent.embeddings, local.spatial_coords))
    rsquared <- c()
    x_coeff <- c()
    y_coeff <- c()
    for (j in 1:ncol(latent.embeddings)){
      s <- summary(lm(calc.df[, j] ~ X + Y, data = calc.df))
      rsquared[j] <- s$r.squared
      x_coeff[j] <- s$coefficients[2, 1]
      y_coeff[j] <- s$coefficients[3, 1]
    }
    best.index <- which.max(rsquared)
    local.coeff <- c(x_coeff[best.index], y_coeff[best.index])
    coeff.list[[i]] <- local.coeff
    assignment.list[i] <- colnames(latent.embeddings)[best.index]
    rsquared.list[i] <- rsquared[best.index]
    full.linear.res <- data.frame(x_coeff = x_coeff, y_coeff = y_coeff, rsquared = rsquared)
    rownames(full.linear.res) <- colnames(latent.embeddings)
    full.res.list[[i]] <- full.linear.res
  }
  coeff.df <- as.data.frame(do.call(rbind, coeff.list))
  colnames(coeff.df) <- c("X", "Y")
  return(list(cell = cell.list, assignment.list = unlist(assignment.list),
              rsquared.list = unlist(rsquared.list), coeff = coeff.df, full.linear = full.res.list
  ))
}

#' main function of the LSGI package
#'
#' @param spatial_coords a dataframe with cell/spot as row and axis (X, Y) as column
#' @param embeddings embeddings used to calculate gradient with
#' @param n.grids.scale number of grids to use equals to (total number of cells)/n.grids.scales
#' @param n.cells.per.meta the size of local group; it is suggested to use 20-50 as a start
#' @return returns local linear gradient information and other basic information of intermediate steps of this analysis
#' @export
#'
local.traj.preprocessing <- function(spatial_coords, embeddings, n.grids.scale = 5, n.cells.per.meta = 25){
  grids <- get.grid.coords(spatial_coords = spatial_coords, n.grids.scale = n.grids.scale)
  dist.to.grid <- spa_vectorized_pdist(A = as.matrix(spatial_coords), B = as.matrix(grids))
  ts <- calc.local.linear(grids = grids, dist.to.grid = dist.to.grid, latent.embeddings = embeddings, spatial_coords = spatial_coords, n.cells.per.meta = n.cells.per.meta)
  grid.info <- cbind(grids, ts$coeff, ts$rsquared.list, ts$assignment.list)
  colnames(grid.info) <- c("X", "Y", "vx", "vy", "R_squared", "Assignment")
  grid.info$Assignment <- factor(x = grid.info$Assignment, levels = colnames(embeddings))
  grid.info <- optimize.arrow(grid.info = grid.info)
  return(list(grid.info = grid.info, local.linear.info = ts, grids = grids,
              dist.to.grid = dist.to.grid, spatial_coords = spatial_coords,
              embeddings = embeddings))
}


#' to optimize arrow length for visualization
#'
#' @param grid.info part of the output of local.traj.preprocessing
#' @param scale.factor scale factor for lengths of arrows; 50 is good for visium data as tested
#' @return returns updated grid.info

optimize.arrow <- function(grid.info, scale.factor = 50){
  grid.info$qsum <- sqrt(grid.info$vx^2 + grid.info$vy^2)
  dim.scale <- mean(max(grid.info$X) - min(grid.info$X),
                    max(grid.info$Y) - min(grid.info$Y))
  grid.info$sf <- (dim.scale/scale.factor)/grid.info$qsum
  grid.info$vx.u <- grid.info$vx * grid.info$sf
  grid.info$vy.u <- grid.info$vy * grid.info$sf
  return(grid.info)
}

#' to gather grid level information from the LSGI result file.
#'
#' @param info LSGI output
#' @return returns a dataframe with all gradient (for every factor) at all locations (for every grid point)

get.ind.rsqrs <- function(info){
  res <- info[["local.linear.info"]][["full.linear"]]
  n.grids <- length(res)
  n.fctrs <- nrow(res[[1]])
  res.df <- do.call(rbind, res)
  colnames(res.df)[c(1, 2)] <- c("vx", "vy")
  res.df$fctr <- rep(colnames(info[["embeddings"]]), n.grids)
  res.df$grid <- rep(paste0("grid_", 1:n.grids), each = n.fctrs)
  res.df$grid <- factor(res.df$grid, levels = unique(res.df$grid))
  res.df$X <- rep(info$grid.info[, "X"], each = n.fctrs)
  res.df$Y <- rep(info$grid.info[, "Y"], each = n.fctrs)
  res.df <- optimize.arrow(res.df)
  return(res.df)
}

#' to plot individual factor superimposed with embeddings/levels
#'
#' @param info LSGI (local.traj.preprocessing) output
#' @param fctr factor to plot
#' @param r_squared_thresh R-squared threshold: only plot gradient with higher R-squared than this value
#' @param arrow.length.scale for adjusting arrow length in a flexible way
#' @return returns a plot with individual factor superimposed with embeddings/levels
#' @export
#'
plt.factor.gradient.ind <- function(info, fctr, r_squared_thresh, arrow.length.scale = 1){

  lin.res.df <- get.ind.rsqrs(info)
  lin.res.df <- na.omit(lin.res.df)
  grid.info <- info$grid.info
  embd <- info[["embeddings"]]
  spatial_coords <- info$spatial_coords
  spatial_coords$factor <- embd[, fctr]
  p <- ggplot(spatial_coords, aes(x=X, y=Y, color = factor)) +
    geom_point(size=2, shape = 20, stroke = 0) +
    scale_color_viridis(direction = -1)+
    geom_segment(data = lin.res.df[lin.res.df$fctr == fctr & lin.res.df$rsquared > r_squared_thresh, ],
                 aes(xend = X + vx.u*arrow.length.scale, yend = Y + vy.u*arrow.length.scale, fill = NULL),
                 color = "darkred",
                 linewidth = 0.4, arrow = arrow(length = unit(0.1, "cm"))) +
    theme_classic() +
    theme(axis.text.x = element_text(face = "bold", color = "black",
                                     size = 12, angle = 0, hjust = 1),
          axis.text.y = element_text(face = "bold", color = "black",
                                     size = 12, angle = 0))
  return(p)
}

#' to plot multiple factors together on the spatial map
#'
#' @param info LSGI (local.traj.preprocessing) output
#' @param sel.factors factor to plot; by default all
#' @param minimum.fctr if some factors have too few gradients, ignore it
#' @param r_squared_thresh R-squared threshold: only plot gradient with higher R-squared than this value
#' @param arrow.length.scale for adjusting arrow length in a flexible way
#' @return returns a plot with multiple selected factors on the spatial map
#' @export
#'
plt.factors.gradient.ind <- function(info, r_squared_thresh, sel.factors = NULL, minimum.fctr = 5,
                                      arrow.length.scale = 1){

  lin.res.df <- get.ind.rsqrs(info)
  lin.res.df <- na.omit(lin.res.df)
  lin.res.df <- lin.res.df[lin.res.df$rsquared > r_squared_thresh, ]
  if (!is.null(sel.factors)){
    lin.res.df <- lin.res.df[lin.res.df$fctr %in% sel.factors, ]
  }
  lin.res.df <- lin.res.df %>%
    group_by(fctr) %>%
    filter(n() >= minimum.fctr) %>%
    ungroup()
  grid.info <- info$grid.info

  spatial_coords <- info$spatial_coords
  p <- ggplot(spatial_coords, aes(x=X, y=Y)) +
    geom_point(size=2, shape = 20, stroke = 0.1, color = "lightgrey") +
    geom_segment(data = lin.res.df, aes(xend = X + vx.u*arrow.length.scale,
                                        yend = Y + vy.u*arrow.length.scale, color = fctr),
                 linewidth = 0.4, arrow = arrow(length = unit(0.1, "cm"))) +
    scale_color_brewer(palette = "Paired") +
    theme_classic() +
    theme(axis.text.x = element_text(face = "bold", color = "black",
                                     size = 12, angle = 0, hjust = 1),
          axis.text.y = element_text(face = "bold", color = "black",
                                     size = 12, angle = 0))
  return(p)
}


#' to plot heatmap of distance between gradated factors
#'
#' @param dist.mat output of avg.dist.calc
#'
#' @return returns a heatmap of distance between gradated factors
#' @export
plt.dist.heat <- function(dist.mat){
  mat2plot <- log2(reshape2::acast(dist.mat, formula = grid.1 ~ grid.2))
  mat2plot <- mat2plot[order(as.numeric(gsub(".*_","",colnames(mat2plot)))),
                       order(as.numeric(gsub(".*_","",colnames(mat2plot))))]
  p1 <- ComplexHeatmap::Heatmap(mat2plot, cluster_rows = F, col = rev(viridis(30)),
                                cluster_columns = F, heatmap_legend_param = list(
                                  title = "log2(distance)"))
  return(p1)
}

#' to calculate distance between gradated factors
#'
#' @param info LSGI (local.traj.preprocessing) output
#' @param minimum.fctr if some factors have too few gradients, ignore it
#' @param r_squared_thresh R-squared threshold: only plot gradient with higher R-squared than this value
#' @return returns a heatmap of distance between gradated factors
#' @export
avg.dist.calc <- function(info, r_squared_thresh = 0.6, minimum.fctr = 10){
  lin.res.df <- get.ind.rsqrs(info)
  lin.res.df <- na.omit(lin.res.df)
  lin.res.df <- lin.res.df[lin.res.df$rsquared > r_squared_thresh, ]
  lin.res.df$rowname.store <- rownames(lin.res.df)
  lin.res.df <- lin.res.df %>%
    group_by(fctr) %>%
    filter(n() >= minimum.fctr) %>%
    ungroup()
  lin.res.df <- as.data.frame(lin.res.df)
  rownames(lin.res.df) <- lin.res.df$rowname.store
  coords.df <- lin.res.df[, c("X", "Y")]
  dist.mat <- spa_vectorized_pdist(A = as.matrix(coords.df),
                                   B = as.matrix(coords.df))
  diag(dist.mat) <- NA
  dist.mat <- reshape2::melt(dist.mat, na.rm=T)
  dist.mat$grid.1 <- lin.res.df[dist.mat$Var1, "fctr"]
  dist.mat$grid.2 <- lin.res.df[dist.mat$Var2, "fctr"]
  dist.mat <- dist.mat[order(as.numeric(gsub(".*_","",dist.mat$grid.1))), ]
  res.list <- list()
  for (fctr in unique(dist.mat$Var1)){
    res <- dist.mat[dist.mat$Var1 == fctr, ] %>%
      group_by(grid.2) %>%
      slice_min(order_by = value, n = 1)
    res.list[[fctr]] <- res
  }
  res.df <- do.call(rbind, res.list)
  for (i in 1:nrow(res.df)){
    a <- res.df[i, "grid.1"]
    b <- res.df[i, "grid.2"]
  }

  res.df$dist.category <- paste0(res.df$grid.1, ".", res.df$grid.2)
  res.df <- res.df %>%
    group_by(dist.category) %>%
    summarise(mean_value = mean(value))
  res.df <- as.data.frame(res.df)

  ret.df <- as.data.frame(do.call(rbind, strsplit(res.df$dist.category, "[.]")))
  colnames(ret.df) <- c("grid.1", 'grid.2')
  ret.df$distance <- res.df$mean_value
  ret.df <- ret.df[order(as.numeric(gsub(".*_","",ret.df$grid.1))), ]
  return(ret.df)
}

