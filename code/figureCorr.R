# ==================================
# = process data for GWAS analysis =
# ==================================

args <- commandArgs(TRUE) # args <- c("../data/OSDGRPmating.csv")
library("RColorBrewer")
library("grid")
library("diagram")

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

cross.data <- read.csv(args[4], header = TRUE, as.is = TRUE)

# plot distribution
# ============================================================

file.width = 178
cairo_pdf(file = args[1], width = file.width/25.4, height = file.width/25.4 * 1/2, family = args[2])
par(las = 1, tcl = -0.2, mai = c(0.05, 0.05, 0.01, 0.01)*file.width/25.4, ps = 7, lwd = 0.5, xpd = F, mfrow = c(2, 4))

# female -> male
# ============================================================

common.line.order <- intersect(intersect(names(os.female.ratio.1hr.line.mean), names(os.male.ratio.30m.line.mean)), names(dgrp.self.ratio.1hr.line.mean))

plot(os.male.ratio.30m.line.mean[common.line.order], os.female.ratio.1hr.line.mean[common.line.order], axes = FALSE, xlim = c(0, 1.05), ylim = c(0, 1.15), xlab = "", ylab = "", pch = 3, cex = 0.5)
abline(a = 0, b = 1, xpd = FALSE)

axis(side = 1, mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "DGRP     x O/S", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(0.52, -0.264, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.80, -0.264, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)

axis(side = 2, mgp = c(1, 0.4, 0), at = seq(0, 1, 0.2), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "O/S     x DGRP", mgp = c(1.6, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
text(-0.27, 0.49, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)
text(-0.27, 0.86, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)

box(bty = "l", lwd = 0.5)

straightarrow(c(0.45, 1.17), c(0.55, 1.17), arr.pos = 1, lty = 1, lwd = 0.5, arr.length = 0.07, arr.width = 0.07, xpd = TRUE)
text(0.40, 1.17, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.62, 1.17, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)

xy.cor <- cor.test(os.male.ratio.30m.line.mean[common.line.order], os.female.ratio.1hr.line.mean[common.line.order])

cor.est <- formatC(xy.cor$estimate, format = "f", digits = 2)
cor.p <- formatC(xy.cor$p.value, format = "f", digits = 2)

text(0.5, 1.07, parse(text = paste("paste(italic(r), \" = \", ", cor.est, ", \" (\", italic(P), \" = \", ", cor.p,  ", \")\" )")), cex = 7/par("ps")/par("cex"))

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# os male vs dgrp
# ============================================================

common.line.order <- intersect(intersect(names(os.female.ratio.1hr.line.mean), names(os.male.ratio.30m.line.mean)), names(dgrp.self.ratio.1hr.line.mean))

plot(os.male.ratio.30m.line.mean[common.line.order], dgrp.self.ratio.1hr.line.mean[common.line.order], axes = FALSE, xlim = c(0, 1.05), ylim = c(0, 1.15), xlab = "", ylab = "", pch = 3, cex = 0.5)
abline(a = 0, b = 1, xpd = FALSE)

axis(side = 1, mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "DGRP     x O/S", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(0.52, -0.264, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.80, -0.264, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)

axis(side = 2, mgp = c(1, 0.4, 0), at = seq(0, 1, 0.2), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "DGRP     x DGRP", mgp = c(1.6, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
text(-0.27, 0.54, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)
text(-0.27, 0.90, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)

box(bty = "l", lwd = 0.5)

straightarrow(c(0.45, 1.17), c(0.55, 1.17), arr.pos = 1, lty = 1, lwd = 0.5, arr.length = 0.07, arr.width = 0.07, xpd = TRUE)
text(0.40, 1.17, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.68, 1.17, "DGRP", cex = 7/par("ps")/par("cex"), xpd = TRUE)

xy.cor <- cor.test(os.male.ratio.30m.line.mean[common.line.order], dgrp.self.ratio.1hr.line.mean[common.line.order])

cor.est <- formatC(xy.cor$estimate, format = "f", digits = 2)
cor.p <- formatC(xy.cor$p.value, format = "e", digits = 2)

text(0.5, 1.07, parse(text = paste("paste(italic(r), \" = \", ", cor.est, ", \" (\", italic(P), \" = \", ", cor.p,  ", \")\" )")), cex = 7/par("ps")/par("cex"), xpd = T)

text(grconvertX(0.05 + file.width/25.4/4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# os female vs dgrp
# ============================================================

common.line.order <- intersect(intersect(names(os.female.ratio.1hr.line.mean), names(os.male.ratio.30m.line.mean)), names(dgrp.self.ratio.1hr.line.mean))

plot(os.female.ratio.1hr.line.mean[common.line.order], dgrp.self.ratio.1hr.line.mean[common.line.order], axes = FALSE, xlim = c(0, 1.05), ylim = c(0, 1.15), xlab = "", ylab = "", pch = 3, cex = 0.5)
abline(a = 0, b = 1, xpd = FALSE)

axis(side = 1, mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "O/S     x DGRP", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(0.45, -0.264, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.80, -0.264, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)

axis(side = 2, mgp = c(1, 0.4, 0), at = seq(0, 1, 0.2), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "DGRP     x DGRP", mgp = c(1.6, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
text(-0.27, 0.54, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)
text(-0.27, 0.90, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)

box(bty = "l", lwd = 0.5)

straightarrow(c(0.45, 1.17), c(0.55, 1.17), arr.pos = 1, lty = 1, lwd = 0.5, arr.length = 0.07, arr.width = 0.07, xpd = TRUE)
text(0.40, 1.17, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.68, 1.17, "DGRP", cex = 7/par("ps")/par("cex"), xpd = TRUE)

xy.cor <- cor.test(os.female.ratio.1hr.line.mean[common.line.order], dgrp.self.ratio.1hr.line.mean[common.line.order])

cor.est <- formatC(xy.cor$estimate, format = "f", digits = 2)
cor.p <- formatC(xy.cor$p.value, format = "e", digits = 2)

text(0.5, 1.07, parse(text = paste("paste(italic(r), \" = \", ", cor.est, ", \" (\", italic(P), \" = \", ", cor.p,  ", \")\" )")), cex = 7/par("ps")/par("cex"), xpd = T)

text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# additive
# ============================================================

common.line.order <- intersect(intersect(names(os.female.ratio.1hr.line.mean), names(os.male.ratio.30m.line.mean)), names(dgrp.self.ratio.1hr.line.mean))

plot((os.female.ratio.1hr.line.mean[common.line.order] + os.male.ratio.30m.line.mean[common.line.order])/2, dgrp.self.ratio.1hr.line.mean[common.line.order], axes = FALSE, xlim = c(0, 1.05), ylim = c(0, 1.15), xlab = "", ylab = "", pch = 3, cex = 0.5)
abline(a = 0, b = 1, xpd = FALSE)

axis(side = 1, mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "(    +    )/2", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(0.43, -0.264, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.55, -0.264, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)

axis(side = 2, mgp = c(1, 0.4, 0), at = seq(0, 1, 0.2), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "DGRP     x DGRP", mgp = c(1.6, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
text(-0.27, 0.54, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)
text(-0.27, 0.90, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)

box(bty = "l", lwd = 0.5)

straightarrow(c(0.45, 1.17), c(0.55, 1.17), arr.pos = 1, lty = 1, lwd = 0.5, arr.length = 0.07, arr.width = 0.07, xpd = TRUE)
text(0.28, 1.17, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.34, 1.17, "+", cex = 7/par("ps")/par("cex"), xpd = TRUE)
text(0.40, 1.17, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.68, 1.17, "DGRP", cex = 7/par("ps")/par("cex"), xpd = TRUE)

xy.cor <- cor.test((os.female.ratio.1hr.line.mean[common.line.order] + os.male.ratio.30m.line.mean[common.line.order])/2, dgrp.self.ratio.1hr.line.mean[common.line.order])

cor.est <- formatC(xy.cor$estimate, format = "f", digits = 2)
cor.p <- formatC(xy.cor$p.value, format = "e", digits = 2)

text(0.5, 1.07, parse(text = paste("paste(italic(r), \" = \", ", cor.est, ", \" (\", italic(P), \" = \", ", cor.p,  ", \")\" )")), cex = 7/par("ps")/par("cex"), xpd = T)

text(grconvertX(0.05 + file.width/25.4/4*3, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("d")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


# multiplicative
# ============================================================

common.line.order <- intersect(intersect(names(os.female.ratio.1hr.line.mean), names(os.male.ratio.30m.line.mean)), names(dgrp.self.ratio.1hr.line.mean))

plot((os.female.ratio.1hr.line.mean[common.line.order] * os.male.ratio.30m.line.mean[common.line.order]), dgrp.self.ratio.1hr.line.mean[common.line.order], axes = FALSE, xlim = c(0, 1.05), ylim = c(0, 1.15), xlab = "", ylab = "", pch = 3, cex = 0.5)
abline(a = 0, b = 1, xpd = FALSE)

axis(side = 1, mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = expression(""%.%""), mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(0.48, -0.264, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.57, -0.264, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)

axis(side = 2, mgp = c(1, 0.4, 0), at = seq(0, 1, 0.2), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "DGRP     x DGRP", mgp = c(1.6, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
text(-0.27, 0.54, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)
text(-0.27, 0.90, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)

box(bty = "l", lwd = 0.5)

straightarrow(c(0.45, 1.17), c(0.55, 1.17), arr.pos = 1, lty = 1, lwd = 0.5, arr.length = 0.07, arr.width = 0.07, xpd = TRUE)
text(0.28, 1.17, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.34, 1.17, expression(""%.%""), cex = 7/par("ps")/par("cex"), xpd = TRUE)
text(0.40, 1.17, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.68, 1.17, "DGRP", cex = 7/par("ps")/par("cex"), xpd = TRUE)

xy.cor <- cor.test((os.female.ratio.1hr.line.mean[common.line.order] * os.male.ratio.30m.line.mean[common.line.order]), dgrp.self.ratio.1hr.line.mean[common.line.order])

cor.est <- formatC(xy.cor$estimate, format = "f", digits = 2)
cor.p <- formatC(xy.cor$p.value, format = "e", digits = 2)

text(0.5, 1.07, parse(text = paste("paste(italic(r), \" = \", ", cor.est, ", \" (\", italic(P), \" = \", ", cor.p,  ", \")\" )")), cex = 7/par("ps")/par("cex"), xpd = T)

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("e")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# multiplicative, log scale
# ============================================================

common.line.order <- intersect(intersect(names(os.female.ratio.1hr.line.mean), names(os.male.ratio.30m.line.mean)), names(dgrp.self.ratio.1hr.line.mean))

plot(log(os.female.ratio.1hr.line.mean[common.line.order] * os.male.ratio.30m.line.mean[common.line.order]), log(dgrp.self.ratio.1hr.line.mean[common.line.order]), axes = FALSE, xlim = c(-4, 0.5), ylim = c(-4, 0.5), xlab = "", ylab = "", pch = 3, cex = 0.5)
abline(a = 0, b = 1, xpd = FALSE)

axis(side = 1, mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "ln    + ln     ", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(-2.1, -5.01, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(-1.2, -5.01, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)

axis(side = 2, mgp = c(1, 0.4, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "ln(DGRP     x DGRP    )", mgp = c(1.6, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
text(-5.15, -1.9, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)
text(-5.15, -0.5, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)

box(bty = "l", lwd = 0.5)

straightarrow(c(-2.15, 0.6), c(-1.7, 0.6), arr.pos = 1, lty = 1, lwd = 0.5, arr.length = 0.07, arr.width = 0.07, xpd = TRUE)
text(-2.9, 0.6, "ln    + ln     ", cex = 7/par("ps")/par("cex"), xpd = TRUE)
text(-3.25, 0.6, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(-2.35, 0.6, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)

text(-1.1, 0.6, "DGRP", cex = 7/par("ps")/par("cex"), xpd = TRUE)

x.var <- log(os.female.ratio.1hr.line.mean[common.line.order] * os.male.ratio.30m.line.mean[common.line.order])
y.var <- log(dgrp.self.ratio.1hr.line.mean[common.line.order])

bad.sample <- which(x.var < log(1e-6) | y.var < log(1e-6))


xy.cor <- cor.test(x.var[-bad.sample], y.var[-bad.sample])

cor.est <- formatC(xy.cor$estimate, format = "f", digits = 2)
cor.p <- formatC(xy.cor$p.value, format = "e", digits = 2)

text(-1.8, 0.2, parse(text = paste("paste(italic(r), \" = \",\"", cor.est, "\", \" (\", italic(P), \" = \", ", cor.p,  ", \")\" )", sep = "")), cex = 7/par("ps")/par("cex"), xpd = T)

text(grconvertX(0.05 + file.width/25.4/4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("f")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


# additive, cross
# ============================================================

plot((cross.data[, 5] + cross.data[, 7])/2, cross.data[, 3], axes = FALSE, xlim = c(0, 1.05), ylim = c(0, 1.15), xlab = "", ylab = "", pch = 3, cex = 0.5)
abline(a = 0, b = 1, xpd = FALSE)

axis(side = 1, mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "(    +    )/2", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(0.43, -0.264, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.55, -0.264, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)

axis(side = 2, mgp = c(1, 0.4, 0), at = seq(0, 1, 0.2), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "DGRP     x DGRP", mgp = c(1.6, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
text(-0.27, 0.54, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)
text(-0.27, 0.90, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)

box(bty = "l", lwd = 0.5)

straightarrow(c(0.45, 1.17), c(0.55, 1.17), arr.pos = 1, lty = 1, lwd = 0.5, arr.length = 0.07, arr.width = 0.07, xpd = TRUE)
text(0.28, 1.17, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.34, 1.17, "+", cex = 7/par("ps")/par("cex"), xpd = TRUE)
text(0.40, 1.17, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.78, 1.17, "DGRP cross", cex = 7/par("ps")/par("cex"), xpd = TRUE)

xy.cor <- cor.test((cross.data[, 5] + cross.data[, 7])/2, cross.data[, 3])

cor.est <- formatC(xy.cor$estimate, format = "f", digits = 2)
cor.p <- formatC(xy.cor$p.value, format = "e", digits = 2)

text(0.5, 1.07, parse(text = paste("paste(italic(r), \" = \", ", cor.est, ", \" (\", italic(P), \" = \", ", cor.p,  ", \")\" )")), cex = 7/par("ps")/par("cex"), xpd = T)

text(grconvertX(0.05 + file.width/25.4/4*2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("g")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# multiplicative, cross
# ============================================================

plot(cross.data[, 5] * cross.data[, 7], cross.data[, 3], axes = FALSE, xlim = c(0, 1.05), ylim = c(0, 1.15), xlab = "", ylab = "", pch = 3, cex = 0.5)
abline(a = 0, b = 1, xpd = FALSE)

axis(side = 1, mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = expression(""%.%""), mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(0.48, -0.264, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.57, -0.264, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)

axis(side = 2, mgp = c(1, 0.4, 0), at = seq(0, 1, 0.2), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = "DGRP     x DGRP", mgp = c(1.6, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
text(-0.27, 0.54, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)
text(-0.27, 0.90, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)

box(bty = "l", lwd = 0.5)

straightarrow(c(0.45, 1.17), c(0.55, 1.17), arr.pos = 1, lty = 1, lwd = 0.5, arr.length = 0.07, arr.width = 0.07, xpd = TRUE)
text(0.28, 1.17, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.34, 1.17, expression(""%.%""), cex = 7/par("ps")/par("cex"), xpd = TRUE)
text(0.40, 1.17, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.78, 1.17, "DGRP cross", cex = 7/par("ps")/par("cex"), xpd = TRUE)

xy.cor <- cor.test(cross.data[, 5] * cross.data[, 7], cross.data[, 3])

cor.est <- formatC(xy.cor$estimate, format = "f", digits = 2)
cor.p <- formatC(xy.cor$p.value, format = "e", digits = 2)

text(0.5, 1.07, parse(text = paste("paste(italic(r), \" = \", ", cor.est, ", \" (\", italic(P), \" = \", ", cor.p,  ", \")\" )")), cex = 7/par("ps")/par("cex"), xpd = T)

text(grconvertX(0.05 + file.width/25.4/4*3, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("h")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)
dev.off()