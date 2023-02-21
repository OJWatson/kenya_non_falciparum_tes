#' Save Figures
#'
#' @param name Name of figure
#' @param fig ggplot or similar figure object
#' @param width Width of plot in inches. Default = 6
#' @param height Height of plot in inches. Default = 6
#' @param plot_dir Plotting directory. Defaults to "analysis/plots"
#'
#' @importFrom grDevices dev.off pdf
save_figs <- function(name,
                      fig,
                      width = 6,
                      height = 6,
                      plot_dir = file.path(here::here(), "analysis/plots")) {

  dir.create(plot_dir, showWarnings = FALSE)
  fig_path <- function(name) {paste0(plot_dir, "/", name)}

  cowplot::save_plot(filename = fig_path(paste0(name,".png")),
                     plot = fig,
                     base_height = height,
                     base_width = width)

  pdf(file = fig_path(paste0(name,".pdf")), width = width, height = height)
  print(fig)
  dev.off()

}

#' Create path relative to root of project
#'
#' @param path Path to be appended to project root
cp_path <- function(path) {

  file.path(here::here(), path)

}

#' Write stargazer model output to file
#'
#' @param x Stargazer model output as text
#' @param path Path for stargazer to be written to
#' @param cls Number of columns in final stargazer table. Default = NULL,
#'   which works it out based on maximum split size
write_stargazer <- function(x, path, cls = NULL) {

  splitup <- vapply(x, strsplit, "\\s{2,}", FUN.VALUE = vector("list", 1))

  if(is.null(cls)) {
    cls <- max(lengths(splitup))
  }
  tbl <- do.call(rbind, splitup[lengths(splitup) == cls])
  rownames(tbl) <- NULL
  colnames(tbl) <- tbl[1,]
  write.csv(tbl[-1,], path, row.names = FALSE)

}

#' Turn gt summary table into a ggplot image
#'
#' @param x gtsummary plot
#' @param ... Other arguments to gt::gtsave
#' @details Ref: https://github.com/MSKCC-Epi-Bio/bstfun/blob/main/R/as_ggplot.R
as_ggplot <- function(x, ...) {
  # checks ---------------------------------------------------------------------
  if (!inherits(x, c("gt_tbl", "gtsummary"))) {
    stop("`x=` must be a 'gt' or 'gtsummary' table", call. = FALSE)
  }


  # convert gtsummary to gt ----------------------------------------------------
  if (inherits(x, "gtsummary")) {
    x <- gtsummary::as_gt(x)
  }

  # save gt as image -----------------------------------------------------------
  path_gt_table_image <- fs::file_temp(ext = "png")
  gt_table_image <- gt::gtsave(x, filename = path_gt_table_image, ...)

  # save image in ggplot -------------------------------------------------------
  table_img <-
    magick::image_read(path_gt_table_image) %>%
    magick::image_ggplot(interpolate = TRUE)

  table_img
}
