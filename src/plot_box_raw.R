plot_box_raw <- function(dat, outcome, relab = NULL){
  
  if(!is.null(relab)){
    colnames(dat) <- gsub(outcome,unname(relab[outcome]),colnames(dat))
    outcome <- unname(relab[outcome])
  }
  
  p <- ggplot(dat) +
    geom_boxplot(aes(x=diet4,y=.data[[outcome]],fill=exposure),
                 alpha=0.6,
                 colour='black',
                 outlier.shape = NA) +
    ggbeeswarm::geom_beeswarm(aes(x=diet4,y=.data[[outcome]],colour=exposure),
                              dodge.width = 0.75,
                              size = 2) + 
    scale_fill_manual(values=c('#EAEAEA','#BBBBBC'),
                      labels = c(
                        "P/D-" = "Healthy",
                        "P/D+" = "Exposed")) +
    scale_colour_manual(values=c('#A9A9A9','#36454F'),
                        labels = c(
                          "P/D-" = "Healthy",
                          "P/D+" = "Exposed")) +
    scale_x_discrete(labels = c(
      "KD- MIF-" = "Ctrl",
      "KD+ MIF-" ="KD",
      "KD- MIF+" = "MIF",
      "KD+ MIF+" = "KD + MIF"
    )) +
    theme_classic() +
    xlab('') +
    theme(legend.title = element_blank())
  return(p)
  
}
