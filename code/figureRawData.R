# ==================================
# = process data for GWAS analysis =
# ==================================

args <- commandArgs(TRUE) # args <- c("../data/OSDGRPmating.csv")
library("RColorBrewer")

# read raw data
# ============================================================
os.data <- read.csv(args[3], header = TRUE, as.is = TRUE, na.strings = c("NA", ""))
os.data <- subset(os.data, !is.na(series) & n.female == 5 & n.male == 10)

# make calculations
os.data <- within(os.data, {female.line <- factor(female.line); male.line <- factor(male.line); ratio.30m <- as.numeric(ratio.30m); ratio.30m2 <- ifelse(is.na(n.mate.30m), 0, n.mate.30m)/5; n.mate.1hr <- as.numeric(n.mate.1hr); ratio.1hr2 <- ifelse(is.na(n.mate.1hr), 0, n.mate.1hr)/5 })

# replace ratio.30m and ratio.1hr in the original with the new calculations
# there were some errors
os.data <- within(os.data, {ratio.30m <- ratio.30m2; ratio.1hr <- ratio.1hr2;})
dgrp.self <- subset(os.data, as.character(female.line) == as.character(male.line) & as.character(female.line) != "OS")
dgrp.self.ratio.1hr.line.mean <- sort(with(dgrp.self, sapply(split(ratio.1hr, male.line), mean)))
dgrp.self.ratio.30m.line.mean <- sort(with(dgrp.self, sapply(split(ratio.30m, male.line), mean)))
os.female <- subset(os.data, as.character(female.line) == "OS" & as.character(male.line) != "OS")
os.female.ratio.1hr.line.mean <- sort(with(os.female, sapply(split(ratio.1hr, male.line), mean)))
os.female.ratio.30m.line.mean <- sort(with(os.female, sapply(split(ratio.30m, male.line), mean)))
os.male <- subset(os.data, as.character(male.line) == "OS" & as.character(female.line) != "OS")
os.male.ratio.30m.line.mean <- sort(with(os.male, sapply(split(ratio.30m, female.line), mean)))
os.male.ratio.1hr.line.mean <- sort(with(os.male, sapply(split(ratio.1hr, female.line), mean)))

# plot distribution
# ============================================================

file.width = 178
cairo_pdf(file = args[1], width = file.width/25.4, height = file.width/25.4 * 3/4, family = args[2])
par(las = 1, tcl = -0.2, mai = c(0.05, 0.05, 0.07, 0.01)*file.width/25.4, ps = 7, lwd = 0.5, xpd = F, mfrow = c(3, 4))

# DGRP female x O/S male
# ============================================================

barplot(table(with(subset(os.data, as.character(female.line) != "OS" & as.character(male.line) == "OS"), cut(ratio.30m, breaks = seq(-0.1, 1.1, 0.2)))), space = 1, axes = FALSE, names.arg = NA, xlab = "", ylab = "", ylim = c(0, 500), col = brewer.pal(9, "Set1")[1])

axis(side = 1, at = seq(2, 12, 2) - 0.5, label = seq(0, 1, 0.2), mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "Mating success after 0.5 hr", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(1, 0.4, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "Count of replicates", mgp = c(1.8, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(main = "DGRP x O/S\ndistribution of replicates", cex.main = 7/par("ps")/par("cex"), col.main = brewer.pal(9, "Set1")[1])

box(bty = "l", lwd = 0.5)
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


barplot(table(with(subset(os.data, as.character(female.line) != "OS" & as.character(male.line) == "OS"), cut(ratio.1hr, breaks = seq(-0.1, 1.1, 0.2)))), space = 1, axes = FALSE, names.arg = NA, xlab = "", ylab = "", ylim = c(0, 500), col = brewer.pal(9, "Set1")[1])

axis(side = 1, at = seq(2, 12, 2) - 0.5, label = seq(0, 1, 0.2), mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "Mating success after 1 hr", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(1, 0.4, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "Count of replicates", mgp = c(1.8, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(main = "DGRP x O/S\ndistribution of replicates", cex.main = 7/par("ps")/par("cex"), col.main = brewer.pal(9, "Set1")[1])

box(bty = "l", lwd = 0.5)
text(grconvertX(0.05 + file.width/25.4/4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


barplot(hist(os.male.ratio.30m.line.mean, breaks = seq(0, 1, 0.05), plot = F)$counts, space = 0, axes = FALSE, names.arg = NA, xlab = "", ylab = "", ylim = c(0, 40), col = brewer.pal(9, "Set1")[1])

axis(side = 1, seq(0, 20, 4), label = seq(0, 1, 0.2), mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "Mating success after 0.5 hr", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(1, 0.4, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "Count of lines", mgp = c(1.2, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(main = "DGRP x O/S\ndistribution of line means", cex.main = 7/par("ps")/par("cex"), col.main = brewer.pal(9, "Set1")[1])

box(bty = "l", lwd = 0.5)
polygon(c(-5.5, -5.5, 21.5, 21.5), c(-13, 55.5, 55.5, -13), lwd = 1, border = brewer.pal(9, "Set1")[1], xpd = TRUE)

text(grconvertX(0.05  + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


barplot(hist(os.male.ratio.1hr.line.mean, breaks = seq(0, 1, 0.05), plot = F)$counts, space = 0, axes = FALSE, names.arg = NA, xlab = "", ylab = "", ylim = c(0, 40), col = brewer.pal(9, "Set1")[1])

axis(side = 1, seq(0, 20, 4), label = seq(0, 1, 0.2), mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "Mating success after 1 hr", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(1, 0.4, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "Count of lines", mgp = c(1.2, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(main = "DGRP x O/S\ndistribution of line means", cex.main = 7/par("ps")/par("cex"), col.main = brewer.pal(9, "Set1")[1])

box(bty = "l", lwd = 0.5)
text(grconvertX(0.05  + file.width/25.4/4*3, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("d")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# O/S female x DGRP male
# ============================================================

barplot(table(with(subset(os.data, as.character(female.line) == "OS" & as.character(male.line) != "OS"), cut(ratio.30m, breaks = seq(-0.1, 1.1, 0.2)))), space = 1, axes = FALSE, names.arg = NA, xlab = "", ylab = "", ylim = c(0, 500), col = brewer.pal(9, "Set1")[2])

axis(side = 1, at = seq(2, 12, 2) - 0.5, label = seq(0, 1, 0.2), mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "Mating success after 0.5 hr", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(1, 0.4, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "Count of replicates", mgp = c(1.8, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(main = "O/S x DGRP\ndistribution of replicates", cex.main = 7/par("ps")/par("cex"), col.main = brewer.pal(9, "Set1")[2])

box(bty = "l", lwd = 0.5)
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("e")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


barplot(table(with(subset(os.data, as.character(female.line) == "OS" & as.character(male.line) != "OS"), cut(ratio.1hr, breaks = seq(-0.1, 1.1, 0.2)))), space = 1, axes = FALSE, names.arg = NA, xlab = "", ylab = "", ylim = c(0, 500), col = brewer.pal(9, "Set1")[2])

axis(side = 1, at = seq(2, 12, 2) - 0.5, label = seq(0, 1, 0.2), mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "Mating success after 1 hr", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(1, 0.4, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "Count of replicates", mgp = c(1.8, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(main = "O/S x DGRP\ndistribution of replicates", cex.main = 7/par("ps")/par("cex"), col.main = brewer.pal(9, "Set1")[2])

box(bty = "l", lwd = 0.5)
text(grconvertX(0.05 + file.width/25.4/4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("f")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


barplot(hist(os.female.ratio.30m.line.mean, breaks = seq(0, 1, 0.05), plot = F)$counts, space = 0, axes = FALSE, names.arg = NA, xlab = "", ylab = "", ylim = c(0, 50), col = brewer.pal(9, "Set1")[2])

axis(side = 1, seq(0, 20, 4), label = seq(0, 1, 0.2), mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "Mating success after 0.5 hr", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(1, 0.4, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "Count of lines", mgp = c(1.2, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(main = "O/S x DGRP\ndistribution of line means", cex.main = 7/par("ps")/par("cex"), col.main = brewer.pal(9, "Set1")[2])

box(bty = "l", lwd = 0.5)
text(grconvertX(0.05  + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("g")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


barplot(hist(os.female.ratio.1hr.line.mean, breaks = seq(0, 1, 0.05), plot = F)$counts, space = 0, axes = FALSE, names.arg = NA, xlab = "", ylab = "", ylim = c(0, 40), col = brewer.pal(9, "Set1")[2])

axis(side = 1, seq(0, 20, 4), label = seq(0, 1, 0.2), mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "Mating success after 1 hr", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(1, 0.4, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "Count of lines", mgp = c(1.2, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(main = "O/S x DGRP\ndistribution of line means", cex.main = 7/par("ps")/par("cex"), col.main = brewer.pal(9, "Set1")[2])

box(bty = "l", lwd = 0.5)
polygon(c(-5.5, -5.5, 21.5, 21.5), c(-13, 55.5, 55.5, -13), lwd = 1, border = brewer.pal(9, "Set1")[2], xpd = TRUE)

text(grconvertX(0.05  + file.width/25.4/4*3, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("h")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


# DGRP female x DGRP male
# ============================================================

barplot(table(with(subset(os.data, as.character(female.line) != "OS" & as.character(male.line) != "OS" & as.character(female.line) == as.character(male.line)), cut(ratio.30m, breaks = seq(-0.1, 1.1, 0.2)))), space = 1, axes = FALSE, names.arg = NA, xlab = "", ylab = "", ylim = c(0, 500), col = brewer.pal(9, "Set1")[4])

axis(side = 1, at = seq(2, 12, 2) - 0.5, label = seq(0, 1, 0.2), mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "Mating success after 0.5 hr", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(1, 0.4, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "Count of replicates", mgp = c(1.8, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(main = "DGRP x DGRP\ndistribution of replicates", cex.main = 7/par("ps")/par("cex"), col.main = brewer.pal(9, "Set1")[4])

box(bty = "l", lwd = 0.5)
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("i")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


barplot(table(with(subset(os.data, as.character(female.line) != "OS" & as.character(male.line) != "OS" & as.character(female.line) == as.character(male.line)), cut(ratio.1hr, breaks = seq(-0.1, 1.1, 0.2)))), space = 1, axes = FALSE, names.arg = NA, xlab = "", ylab = "", ylim = c(0, 500), col = brewer.pal(9, "Set1")[4])

axis(side = 1, at = seq(2, 12, 2) - 0.5, label = seq(0, 1, 0.2), mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "Mating success after 1 hr", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(1, 0.4, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "Count of replicates", mgp = c(1.8, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(main = "DGRP x DGRP\ndistribution of replicates", cex.main = 7/par("ps")/par("cex"), col.main = brewer.pal(9, "Set1")[4])

box(bty = "l", lwd = 0.5)
text(grconvertX(0.05 + file.width/25.4/4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("j")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


barplot(hist(dgrp.self.ratio.30m.line.mean, breaks = seq(0, 1, 0.05), plot = F)$counts, space = 0, axes = FALSE, names.arg = NA, xlab = "", ylab = "", ylim = c(0, 50), col = brewer.pal(9, "Set1")[4])

axis(side = 1, seq(0, 20, 4), label = seq(0, 1, 0.2), mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "Mating success after 0.5 hr", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(1, 0.4, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "Count of lines", mgp = c(1.2, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(main = "DGRP x DGRP\ndistribution of line means", cex.main = 7/par("ps")/par("cex"), col.main = brewer.pal(9, "Set1")[4])

box(bty = "l", lwd = 0.5)
text(grconvertX(0.05  + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("k")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


barplot(hist(dgrp.self.ratio.1hr.line.mean, breaks = seq(0, 1, 0.05), plot = F)$counts, space = 0, axes = FALSE, names.arg = NA, xlab = "", ylab = "", ylim = c(0, 40), col = brewer.pal(9, "Set1")[4])

axis(side = 1, seq(0, 20, 4), label = seq(0, 1, 0.2), mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "Mating success after 1 hr", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(1, 0.4, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "Count of lines", mgp = c(1.2, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(main = "DGRP x DGRP\ndistribution of line means", cex.main = 7/par("ps")/par("cex"), col.main = brewer.pal(9, "Set1")[4])

box(bty = "l", lwd = 0.5)
polygon(c(-5.5, -5.5, 21.5, 21.5), c(-13, 55.5, 55.5, -13), lwd = 1, border = brewer.pal(9, "Set1")[4], xpd = TRUE)
text(grconvertX(0.05  + file.width/25.4/4*3, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("l")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


dev.off()