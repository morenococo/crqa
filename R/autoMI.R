autoMI = function(x, nbins, maxlag){
  require(gplots)
  
  ami = vector('numeric', length = maxlag+1);
  
  # Compute a histogram of x
  help(hist)
  
  ## reconstruct the sequence of breaks
  seqb = seq(min(x), max(x), length.out = nbins+1);
  
  # extract the counts
  p_i = hist(x, breaks = seqb, plot = F)$counts;
  # p_i = histc(x, edges =  seqb), plot = F)$counts;
  # Normalize to make it a probability
  p_i = p_i / sum(p_i);
  
  # Save the original input vector
  xorig = x;
  
  # Loop over time lags from 0 to maxlag
  tau = 1;
  for (tau in 0:maxlag){
    
    # Make a time delayed version of x
    xtau = circshift(xorig,-tau);
    
    # Remove the last tau elements of this vector and of x
    xtau = xtau[1:(length(xtau)-tau)];
    x    = xorig[1:(length(xorig)-tau)];
    
    # Compute the two-dimensional histogram of x and xtau
    p_ij = hist2d(x, xtau, nbins, show = F)$counts;
    p_ij = p_ij/sum(p_ij)
    
    # Normalize the rows to get the joint probability the x is
    # in bin i and xtau is in bin j.
    # row_sum = sum(p_ij, 2);
    # If any rows sum to zero, replace the sum by 1
    # row_sum(row_sum == 0) = 1;
    # Divide each element with the row sum. 
    # p_ij = diag(1 ./ row_sum) * p_ij;
    # Compute the auto mutual information by summing over all
    # combinations of i and j
    i = 1; j = 1
    for (i in 1:nbins){
      for (j in 1:nbins){
        
        if (p_i[i] != 0 & p_i[j] != 0 & p_ij[i,j] != 0){
          ami[tau+1] = ami[tau+1] + p_ij[i,j] * log(p_ij[i,j] / (p_i[i] * p_i[j]))};
      }
    }
  }
  return(ami)
}
