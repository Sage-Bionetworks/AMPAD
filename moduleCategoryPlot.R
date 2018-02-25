library(plyr)
library(ggplot2)

setwd("~/work/")

abc <- read.csv(file = "Job-53923557810868363728834403.csv")
# Throw in a bonferroni correction
abc$p2 <- abc$fisherPval * nrow(abc)


# Summarize modules by number of "significant" associations
byModule <- aggregate(p2 ~ ModuleNameFull, data = abc, FUN = function(x) sum(x < .05))
byModule <- arrange(byModule, p2)
byModule$ModuleNameFull <- factor(byModule$ModuleNameFull, levels = byModule$ModuleNameFull, ordered = TRUE)

testPlot <- ggplot(byModule) + theme_bw() +
  geom_point(aes(x = p2, y = ModuleNameFull)) +
  xlab("Number of Categories with Statistically Significant Overlap") +
  ylab("Module") +
  ggtitle("Overlaps Between Genes in Modules and Categories\nMultiple Testing Correction, Fisher's Exact Test")
testPlot
ggsave(file = "testPlot.png", testPlot, width = 8, height = 6)



# Summarize top associations by Module and Category. Grab significant associations, sort by odds ratio
def <- ddply(abc, .(ModuleNameFull), function(x) {
  x <- x[which(x$p2 < .05), ]
  if(!nrow(x)) {
    return(NULL)
  }
  x <- arrange(x, fisherOR)
  # x <- x[1:min(5, nrow(x)), ]
  x
})
def <- arrange(def, p2)
ghi <- def[1:75, ]
jkl <- aggregate(category ~ ModuleNameFull, data = ghi, FUN = length)
jkl <- arrange(jkl, -category)
ghi$ModuleNameFull <- factor(ghi$ModuleNameFull, levels = jkl$ModuleNameFull, ordered = TRUE)
mno <- aggregate(ModuleNameFull ~ category, data = ghi, FUN = length)
mno <- arrange(mno, ModuleNameFull)
ghi$category <- factor(ghi$category, levels = mno$category, ordered = TRUE)

testPlot2 <- ggplot(ghi) + theme_bw() +
  geom_bin2d(aes(x = ModuleNameFull, y = category)) +
  ggtitle("Overlaps Between Genes in Modules and Categories\nTop 75 Associations by Odds Ratio") +
  xlab("Module") + ylab("Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")
testPlot2
ggsave(file = "testPlot2.png", testPlot2, width = 12, height = 8)



# Summarize a single module
modToSum <- "aggregateTCXblueTCX"
d1 <- abc[abc$ModuleNameFull == modToSum & abc$p2 < .05, ]
d1 <- arrange(d1, desc(fisherOR))
d1$category <- factor(d1$category, levels = d1$category, ordered = TRUE)
testPlot3 <- ggplot(d1[1:50, ]) + theme_bw() +
  geom_point(aes(x = fisherOR, y = category)) +
  xlab("Fisher Odds Ratio") + ylab("Category") +
  ggtitle(sprintf("Overlaps between Genes in %s and Categories", modToSum, "\nTop 50 Associations"))
testPlot3
ggsave(file = "testPlot3.png", testPlot3, width = 12, height = 8)


