multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

###

path.cafexp <- "~/Documents/IU/hahn_lab_stuff/CAFExp/"

## max.001 <- read.table(paste0(path.cafexp, "results_0.001_p10_max.txt"), header=TRUE)
max.001 <- read.table(paste0(path.cafexp, "results_0.001_unif_max.txt"), header=TRUE)
## max.005 <- read.table(paste0(path.cafexp, "results_0.005_p10_max.txt"), header=TRUE)
max.005 <- read.table(paste0(path.cafexp, "results_0.005_unif_max.txt"), header=TRUE)
## max.01 <- read.table(paste0(path.cafexp, "results_0.01_p10_max.txt"), header=TRUE)
max.01 <- read.table(paste0(path.cafexp, "results_0.01_unif_max.txt"), header=TRUE)
names(max.001)[2] <- "nlnL"
names(max.005)[2] <- "nlnL"
names(max.01)[2] <- "nlnL"

max.001.plot <- ggplot(max.001, aes(x=l, y=nlnL, colour=as.factor(eqfreq))) +
  geom_point() + geom_line() + xlab(expression(lambda)) + ylab("-lnL") + geom_vline(xintercept=0.001, colour="red") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(color="black"),
    axis.text.x = element_text(color="black", size=12),
    axis.text.y = element_text(color="black", size=12),
    legend.key = element_blank()
  ) + labs(colour="Root equibilibrium\ndistribution") + scale_colour_brewer(palette="Dark2")
max.001.plot

max.005.plot <- ggplot(max.005, aes(x=l, y=nlnL, colour=as.factor(eqfreq))) +
  geom_point() + geom_line() + xlab(expression(lambda)) + ylab("-lnL") + geom_vline(xintercept=0.005, colour="red") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(color="black"),
    axis.text.x = element_text(color="black", size=12),
    axis.text.y = element_text(color="black", size=12),
    legend.key = element_blank()
  ) + labs(colour="Root equibilibrium\ndistribution") + scale_colour_brewer(palette="Dark2")
max.005.plot

max.01.plot <- ggplot(max.01, aes(x=l, y=nlnL, colour=as.factor(eqfreq))) +
  geom_point() + geom_line() + xlab(expression(lambda)) + ylab("-lnL") + geom_vline(xintercept=0.01, colour="red") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(color="black"),
    axis.text.x = element_text(color="black", size=12),
    axis.text.y = element_text(color="black", size=12),
    legend.key = element_blank()
  ) + labs(colour="Root equibilibrium\ndistribution") + scale_colour_brewer(palette="Dark2")
max.01.plot

## pdf(paste0(path.cafexp, "lk_surf_0.001_0.005_0.01_p10_max.pdf"), width=10, height=7.5)
pdf(paste0(path.cafexp, "lk_surf_0.001_0.005_0.01_unif_max.pdf"), width=10, height=7.5)
multiplot(max.001.plot, max.005.plot, max.01.plot, cols=1)
dev.off()

## pdf(paste0(path.cafexp, "lk_surf_0.001_max.pdf"), width=10, height=2.5)
## max.001.plot
## dev.off()

## pdf(paste0(path.cafexp, "lk_surf_0.005_max.pdf"), width=10, height=2.5)
## max.005.plot
## dev.off()

## pdf(paste0(path.cafexp, "lk_surf_0.01_max.pdf"), width=10, height=2.5)
## max.01.plot
## dev.off()
