# Sample of ways to play with the assembly summaries

# Note that you should run this script piecemeal - the different blocks of plots replace each other if you run
# the script from the beginning and you'll only see the last one.

# things to try
#  - only look at specific types of errors (e.g., ignore coverage type errors)
#    do the results change?
#  - Can you see things better in log space?

setwd("./")

# Reading the summary file
#col_nums = c(2,4,6,8) # num errors
#col_nums = c(3,5,7,9) # tot error bps
#col_nums = c(3,5) # low/high cov bps
#col_nums = c(2,4) # num errors
#col_nums = c(6,8) # num breakpoint
#col_nums = c(7,9) # tot breakpoint bps
#col_nums = c(2,4,8) # num errors cov, breakpoints
#col_nums = c(3,5,9)
#col_nums = c(2,4,6,8) # num errors
#col_nums = c(3,5,7,9) # tot error bps
#col_nums = 9 # num errors
#col_nums = c(3,5)
#col_nums = c(7)

#col_list = list(c(2,4,6,8), c(3,5,7,9), c(2,4), c(3,5), 6, 7, 8, 9)
col_list = list(c(3,5,7,9), c(4,6,8,10), 3, 4, 5, 6, 7, 8, 9, 10)
titles = list("All Errors", "All Errors (BPs)", "Low Coverage Errors", "Low Coverage Errors (BPs)", "High Coverage Errors", "High Coverage Errors (BPs)", "Insert Size Errors", "Insert Size Errors (BPs)", "Breakpoint Errors", "Breakpoint Errors (BPs)")

args <- commandArgs(TRUE)
assemblies <- strsplit(args[1], ",")
if (length(args) > 1) {
  names <- strsplit(args[2], ",")
} else {
  names <- assemblies
}

length(args)
if (length(args) > 2) {
  'testA'
  output_filename <- paste(args[3], '.pdf', sep='')
} else {
  'testB'
  output_filename <- 'test.pdf'
}
#tables <- list()

pdf(file=output_filename)
par(mfrow=c(1,1))

# Keep track of which descriptions we should use for each plot.
description_count = 1

for (col_nums in col_list) {
  tables = list()
  print(col_nums)
  xrange <- 10
  yrange <- 10
  for (assembly in assemblies[[1]]) {
    summary = read.csv(assembly, head=T, row.names=1, sep="\t", check.names=F)
    errors = summary[,col_nums]

    if (length(col_nums) > 1) {
      totErrors <- rowSums(errors)
    } else {
      totErrors <- errors
    }
    # weight errors by abundance
    #totErrors <- totErrors * summary[,2]

    # original
    contigSizes = summary[order(summary[, 1], decreasing=T), 1] * summary[order(summary[, 1] * summary[, 2], decreasing=T), 2] / summary[order(summary[, 1] * summary[, 2], decreasing=T), 2]
    contigErrors = totErrors[order(summary[,1], decreasing=T)]

    # abundance-based error sorting
    #contigSizes = summary[order(summary[, 1] * summary[, 2], decreasing=T), 1] * summary[order(summary[, 1] * summary[, 2], decreasing=T), 2]
    #contigErrors = totErrors[order(summary[,1] * summary[, 2], decreasing=T)]

    contigSizes = rbind(0, contigSizes)
    contigErrors = rbind(0, contigErrors)

    contigSums = cumsum(as.numeric(contigSizes))
    errorSums = cumsum(as.numeric(contigErrors))

    contigNum = c(rep(1, length(contigSizes)))
    contigNum = cumsum(as.numeric(contigNum))

    assembly_stats = list(summary, errors, totErrors, contigSizes, contigErrors, contigSums, errorSums, contigNum)
    tables <- c(tables, assembly_stats)

    xrange <- max(xrange, errorSums, na.rm=T)
    yrange <- max(yrange, contigSums, na.rm=T)
  }

  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  count = 0
  #pdf(file=output_filename)
  plot(0,0, pch="",
       xlab="# errors",
       ylab="Cumulative assembly",
       ylim=c(0,yrange),
       xlim=c(0,xrange),
       main = titles[[description_count]]
       #col=cbbPalette[[count+1]]
    )

  legend("bottomright", names[[1]], lwd=c(2.5,2.5), col=cbbPalette, box.col=NA)

  for (assembly in assemblies[[1]]) {
    #plot(tables[[2]][[7]], tables[[2]][[6]], #/contigSums.hmp[length(contigSums.hmp)],

    lines(tables[[count * 8 + 7]], tables[[count * 8 + 6]],
          col=cbbPalette[[count+1]])


    #plot(tables[[count * 8 + 7]], tables[[count * 8 + 6]],
    #     type='l', # log="y",
    #     main = paste("FRC", sample, "- Contigs added by length"),
    #     xlab="# errors",
    #     ylab="Cumulative assembly",
    #     ylim=c(0,yrange),
    #     xlim=c(0,xrange),
    #     col=cbbPalette[[count+1]],

         #legend("topleft", assemblies, lwd=c(2.5,2.5), col=cbbPalette, box.col=NA)
         #log='y'
    #)

    count = count + 1
    #par(new=T)
  }

  description_count = description_count + 1
}
#par(new=F)
dev.off()
warnings()
