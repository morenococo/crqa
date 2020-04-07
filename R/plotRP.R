## plot the RP
## need to be improved, need to take out the for() loop for it
## RP  : the RP plot
## par : a list of parameters for the plotting
##     unit:   the unit to mark ticks on the axes
##     labelx: the label of the x-axis  
##     labely: the label of the y-axis  
##     cols  : the color for the recurrent recurrent points
##     pcex   : the size of the dots

.packageName <- 'crqa'

plotRP <- function(RP, 
                   labelx = "Time", labely = "Time", 
                   labelmain = "Recurrence Plot",
                   point_col = "black", pcex = .3, pch = 1,
                   unit = 10, show_ticks = FALSE){
  
  # find the size of the RP
  xdim   = nrow(RP)
  ydim   = ncol(RP)
  
  # transform the RP into a square matrix
  RP = matrix(as.numeric(RP), nrow = xdim, ncol = ydim)
  
  # figure out where the recurrent points are
  ind = which(RP == 1, arr.ind = T)
  
  # create all time series of all possible samples on axis
  tstamp = seq(0, xdim, 1)
  
  # create the shell of the plot
  par(mar = c(3.8, 3.8, 2,2), font.axis = 2, cex.axis = 1,
      font.lab = 2, cex.lab = 1.2)
  plot(tstamp, tstamp, type = "n", 
       xlab = "", ylab = "", main=labelmain,
       xaxt = "n", yaxt = "n",
       font = 2)
  
  # add recurrent points to the plot
  matpoints(ind[,1], ind[,2],  cex = pcex, col = cols, pch = pch) 
  
  # add x- and y-axis labels to the plot (so they're not too far away)
  mtext(labelx, at = mean(tstamp), side = 1, line = 1.1, cex = 1.2, font = 2)
  mtext(labely, at = mean(tstamp), side = 2, line = 1.1, cex = 1.2, font = 2)
  
  # if the user would like x- and y-tickmarks, show and label them
  if (show_ticks == TRUE){
    labax = seq(0, nrow(RP), unit)
    labay = seq(0, nrow(RP), unit)
    las = 0
    mtext(labax, at = labax, side = 1, line = .1, cex = .8, font = 2, las = las)
    mtext(labay, at = labay, side = 2, line = .1, cex = .8, font = 2, las = las)
  }
  
}
