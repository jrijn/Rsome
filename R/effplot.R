#' Plot concentration curve and linear regression.
#'
#' This function plots the untransformed Cq values for each target against the log2([cDNA]).
#' Then it performs linear regression and plots the R^2 and y ~ x formula.
#' 
#' @param tablename: a relex object with an added column of log2([cDNA]) called logConcentration
#' @param groupColumn: defaults to "Target". Defines the grouping column which will define coloring.
#' @param logConcentration: defaults to "logConcentration". Defines the column containing log2([cDNA]) values.
#' @param meanCq: defaults to "meanCq". Defines the column containing meanCq values.
#' @param plotWidth: width of output pdf file
#' @param plotHeight: height of output pfd file
#' @keywords qPCR, relative expression, plot
#' @export eff.plot
#' @examples
#' 
#' cq <- Rsome::cqimport(cqfile)
#' mc <- Rsome::mcimport(
#'   cqimport = cq, 
#'   meltderivative = meltfile)
#'   
#' re <- Rsome::relex(
#'   cq, 
#'   household = "Gapdh",
#'   SDcutoff = 1,
#'   Cqcutoff = 35)
#' 
#' library(tidyr)
#' library(dplyr)
#' 
#' eff <- separate(re, "Sample", c("Exp", "Condition", "Concentration"), sep = "_")
#' eff$Concentration <- as.numeric(eff$Concentration)
#' eff <- eff %>%
#'   filter(Concentration > 0) %>%
#'   mutate(logConcentration = log2(Concentration))
#' 
#' p1 <- eff.plot(
#'   eff,
#'   logConcentration = "logConcentration",
#'   groupColumn = "Target",
#'   meanCq = "meanCq")
#'   

eff.plot <- function(tablename, 
                    logConcentration = "logConcentration", 
                    groupColumn = "Target", 
                    meanCq = "meanCq",
                    plotWidth = 10, 
                    plotHeight = 10){
  
  #Load required packages
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  #Modify stat_smooth function to include formula and R^2 in plot
  
  stat_smooth_func <- function(mapping = NULL, data = NULL,
                               geom = "smooth", position = "identity",
                               ...,
                               method = "auto",
                               formula = y ~ x,
                               se = TRUE,
                               n = 80,
                               span = 0.75,
                               fullrange = FALSE,
                               level = 0.95,
                               method.args = list(),
                               na.rm = FALSE,
                               show.legend = NA,
                               inherit.aes = TRUE,
                               xpos = NULL,
                               ypos = NULL) {
    layer(
      data = data,
      mapping = mapping,
      stat = StatSmoothFunc,
      geom = geom,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
        method = method,
        formula = formula,
        se = se,
        n = n,
        fullrange = fullrange,
        level = level,
        na.rm = na.rm,
        method.args = method.args,
        span = span,
        xpos = xpos,
        ypos = ypos,
        ...
      )
    )
  }
  
  
  StatSmoothFunc <- ggproto("StatSmooth", Stat,
                            
                            setup_params = function(data, params) {
                              # Figure out what type of smoothing to do: loess for small datasets,
                              # gam with a cubic regression basis for large data
                              # This is based on the size of the _largest_ group.
                              if (identical(params$method, "auto")) {
                                max_group <- max(table(data$group))
                                
                                if (max_group < 1000) {
                                  params$method <- "loess"
                                } else {
                                  params$method <- "gam"
                                  params$formula <- y ~ s(x, bs = "cs")
                                }
                              }
                              if (identical(params$method, "gam")) {
                                params$method <- mgcv::gam
                              }
                              
                              params
                            },
                            
                            compute_group = function(data, scales, method = "auto", formula = y~x,
                                                     se = TRUE, n = 80, span = 0.75, fullrange = FALSE,
                                                     xseq = NULL, level = 0.95, method.args = list(),
                                                     na.rm = FALSE, xpos=NULL, ypos=NULL) {
                              if (length(unique(data$x)) < 2) {
                                # Not enough data to perform fit
                                return(data.frame())
                              }
                              
                              if (is.null(data$weight)) data$weight <- 1
                              
                              if (is.null(xseq)) {
                                if (is.integer(data$x)) {
                                  if (fullrange) {
                                    xseq <- scales$x$dimension()
                                  } else {
                                    xseq <- sort(unique(data$x))
                                  }
                                } else {
                                  if (fullrange) {
                                    range <- scales$x$dimension()
                                  } else {
                                    range <- range(data$x, na.rm = TRUE)
                                  }
                                  xseq <- seq(range[1], range[2], length.out = n)
                                }
                              }
                              # Special case span because it's the most commonly used model argument
                              if (identical(method, "loess")) {
                                method.args$span <- span
                              }
                              
                              if (is.character(method)) method <- match.fun(method)
                              
                              base.args <- list(quote(formula), data = quote(data), weights = quote(weight))
                              model <- do.call(method, c(base.args, method.args))
                              
                              m = model
                              eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                                               list(a = format(coef(m)[[1]], digits = 3), 
                                                    b = format(coef(m)[[2]], digits = 3), 
                                                    r2 = format(summary(m)$r.squared, digits = 3)))
                              func_string = as.character(as.expression(eq))
                              
                              if(is.null(xpos)) xpos = min(data$x)*0.9
                              if(is.null(ypos)) ypos = max(data$y)*0.9
                              data.frame(x=xpos, y=ypos, label=func_string)
                              
                            },
                            
                            required_aes = c("x", "y")
  )
  
  #Then plot the Cq values, faceting is set to the target gene.
  
  p1 <- ggplot(eff, aes(x=logConcentration, y=meanCq), group = groupColumn)+
    scale_y_continuous(limits=c(12.5,NA))+
    geom_point()+
    facet_wrap(
      ~ Target,
      scales = "free")+
    stat_smooth(method = "lm", se = FALSE)+
    stat_smooth_func(
      geom = "text", 
      method = "lm", 
      hjust = 0,
      parse = TRUE)+
    theme_classic()
  
  return(p1)
}
