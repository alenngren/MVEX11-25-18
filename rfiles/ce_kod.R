setwd("/Users/alexander/Chalmers/MVEX11-25-18/fig/clark")

{
  setEPS()
  postscript("CE-thomas.eps", width = 6, height = 3.7)
  
  # Set margins: bottom, left, top, right
  par(mar = c(5, 4, 1, 2))  # reduced top margin to remove empty space
  
  ce_bar <- barplot(height = df$score, names.arg = df$type,
                    col = my_colors,
                    xlab = "",
                    ylab = "R värde",
                    main = "",     # No title
                    cex.main = 1.5,
                    ylim = c(0, 1),
                    cex.lab = 1.2,
                    cex.axis = 1.2
  )
  
  abline(h = ce_x$statistic, lwd = 3, lty = 2, col = "black")
  abline(h = 1, lwd = 3, lty = 2, col = "red")
  
  legend("topleft", legend = c("Data", "CSR-line"), 
         lty = 2, lwd = 3, col = c("black", "red"), bty = "n",
         cex = 1.2)
  
  text(ce_bar, df$score,
       labels = paste("R=", round(df$score, 2), sep = ""),
       pos = 1, cex = 1.2, col = "black")
  
  dev.off()
}

# ce_inhom <- clarkevans.test(X_inhom, alternative="clustered") 
# ce_intlgcp <- clarkevans.test(Xlgcp, alternative="clustered") # intensitetstyrd Thomas
ce_hier <- clarkevans.test(pp_hier, alternative="clustered") 
ce_intp <- clarkevans.test(pp_intp, alternative="clustered") # intensitetstyrd Thomas


my_colors_proc <- c( jet_colors[2], jet_colors[13], jet_colors[10]) # jet_colors[6],
df_proc <- data.frame( # utan ves13 
  score = c(ce_iter10$statistic, ce_intp$statistic, ce_hier$statistic),
  type = c('Iterativ (10)', 'Hierarkisk', 'Intensitetstyrd'),
  stringsAsFactors = FALSE # important to avoid factor issues!
)

{
  setEPS()
  postscript("CE-processes1.eps", width = 6, height = 3.7) # square , width = 6, height = 3.7
  par(mar = c(5, 4, 1, 2))  # reduced top margin to remove empty space
  # Create the barplot with individual bar colors
  ce_bar <- barplot(height=df_proc$score, names.arg=df_proc$type,
                    col=my_colors_proc,
                    xlab="",
                    ylab="R värde",
                    main= "" , # Clark-Evans test för punktprocesser",
                    cex.main=1.5,  # Increase this value for a larger title
                    ylim=c(0,1)
  )
  
  abline(h = ce_x$statistic, lwd=3, lty=2, col="black")
  abline(h = 1, lwd=3, lty=2, col="red")
  legend("topleft", legend = c("Data", "CSR-line"), 
         lty = 2, lwd = 3, col = c("black","red"), bty = "n",
         cex = 1.2)
  
  # Add text labels for R values
  text(ce_bar, df_proc$score,
       labels=paste("R=", round(df_proc$score, 2), sep=""),
       pos=1, cex=1.2, col=c("white", "white", "black"))
  
  dev.off()
}


dirichlet(X)
plot(dirichlet(X), main='')
