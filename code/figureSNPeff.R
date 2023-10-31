# ===========================
# = figure for GWAS effects =
# ===========================

args <- commandArgs(TRUE) # args <- c("../data/DGRPxOSeff.csv")
library("RColorBrewer")
library("grid")
library("diagram")

# read raw data
# ============================================================
eff <- read.csv(args[3], header = T, as.is = T)

# plot effects
# ============================================================

file.width = 89*2
cairo_pdf(file = args[1], width = file.width/25.4, height = file.width/25.4/4, family = args[2])
par(las = 1, tcl = -0.2, mai = c(0.045, 0.07*3/4, 0.07*3/4, 0.015)*file.width/25.4, ps = 7, lwd = 0.5, xpd = F, mfrow = c(1, 4))

# effect between female and male
# ============================================================

plot(eff$OSmale.Eff, eff$OSfemale.Eff, axes = FALSE, xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2), xlab = "", ylab = "", pch = 3, cex = 0.5)
abline(a = 0, b = 1, xpd = FALSE)

axis(side = 1, mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "DGRP     x O/S     SNP effect", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(-0.09, -0.307, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.02, -0.307, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)

axis(side = 2, mgp = c(1, 0.4, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(ylab = "O/S    x DGRP     SNP effect", mgp = c(1.6, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
text(-0.307, -0.135, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)
text(-0.307, 0.02, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)

box(bty = "l", lwd = 0.5)

straightarrow(c(-0.02, 0.32), c(0.02, 0.32), arr.pos = 1, lty = 1, lwd = 0.5, arr.length = 0.07, arr.width = 0.07, xpd = TRUE)
segments(-0.02, 0.32, 0.02, 0.32, xpd = TRUE)
text(-0.04, 0.32, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.05, 0.32, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)

xy.cor <- cor.test(eff$OSmale.Eff, eff$OSfemale.Eff)
cor.est <- formatC(xy.cor$estimate, format = "f", digits = 2)
cor.p <- formatC(xy.cor$p.value, format = "f", digits = 2)
print(cor.p)

text(0, 0.28, parse(text = paste("paste(italic(r), \" = \", ", cor.est, ", \" (\", italic(P), \" = \", ", cor.p,  ", \")\" )")), cex = 7/par("ps")/par("cex"), xpd = TRUE)

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# effect between female and dgrp
# ============================================================

plot(eff$OSmale.Eff, eff$DGRP.Eff, axes = FALSE, xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2), xlab = "", ylab = "", pch = 3, cex = 0.5)
abline(a = 0, b = 1, xpd = FALSE)

axis(side = 1, mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "DGRP     x O/S     SNP effect", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(-0.09, -0.307, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.02, -0.307, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)

axis(side = 2, mgp = c(1, 0.4, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(ylab = "DGRP     x DGRP     SNP effect", mgp = c(1.6, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
text(-0.307, -0.121, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)
text(-0.307, 0.04, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)

box(bty = "l", lwd = 0.5)

straightarrow(c(-0.02, 0.32), c(0.02, 0.32), arr.pos = 1, lty = 1, lwd = 0.5, arr.length = 0.07, arr.width = 0.07, xpd = TRUE)
segments(-0.02, 0.32, 0.02, 0.32, xpd = TRUE)
text(-0.04, 0.32, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.09, 0.32, "DGRP", cex = 7/par("ps")/par("cex"), xpd = TRUE)

xy.cor <- cor.test(eff$OSmale.Eff, eff$DGRP.Eff)
cor.est <- formatC(xy.cor$estimate, format = "f", digits = 2)
cor.p <- formatC(xy.cor$p.value, format = "e", digits = 2)
print(cor.p)

text(0, 0.28, parse(text = paste("paste(italic(r), \" = \", ", cor.est, ", \" (\", italic(P), \" = \", ", cor.p,  ", \")\" )")), cex = 7/par("ps")/par("cex"), xpd = TRUE)

text(grconvertX(0.05 + file.width/25.4/4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# effect between male and dgrp
# ============================================================

plot(eff$OSfemale.Eff, eff$DGRP.Eff, axes = FALSE, xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2), xlab = "", ylab = "", pch = 3, cex = 0.5)
abline(a = 0, b = 1, xpd = FALSE)

axis(side = 1, mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "O/S     x DGRP     SNP effect", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(-0.115, -0.307, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.02, -0.307, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)

axis(side = 2, mgp = c(1, 0.4, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(ylab = "DGRP     x DGRP     SNP effect", mgp = c(1.6, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
text(-0.307, -0.121, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)
text(-0.307, 0.04, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)

box(bty = "l", lwd = 0.5)

straightarrow(c(-0.02, 0.32), c(0.02, 0.32), arr.pos = 1, lty = 1, lwd = 0.5, arr.length = 0.07, arr.width = 0.07, xpd = TRUE)
segments(-0.02, 0.32, 0.02, 0.32, xpd = TRUE)
text(-0.04, 0.32, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.09, 0.32, "DGRP", cex = 7/par("ps")/par("cex"), xpd = TRUE)

xy.cor <- cor.test(eff$OSfemale.Eff, eff$DGRP.Eff)
cor.est <- formatC(xy.cor$estimate, format = "f", digits = 2)
cor.p <- formatC(xy.cor$p.value, format = "e", digits = 2)
print(cor.p)

text(0, 0.28, parse(text = paste("paste(italic(r), \" = \", ", cor.est, ", \" (\", italic(P), \" = \", ", cor.p,  ", \")\" )")), cex = 7/par("ps")/par("cex"), xpd = TRUE)

text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# average
# ============================================================

plot((eff$OSfemale.Eff + eff$OSmale.Eff)/2, eff$DGRP.Eff, axes = FALSE, xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2), xlab = "", ylab = "", pch = 3, cex = 0.5)
abline(a = 0, b = 1, xpd = FALSE)

axis(side = 1, mgp = c(0.8, 0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = "SNP effect (    +    )/2", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(0.035, -0.307, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.085, -0.307, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)

axis(side = 2, mgp = c(1, 0.4, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(ylab = "DGRP     x DGRP     SNP effect", mgp = c(1.6, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
text(-0.307, -0.121, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)
text(-0.307, 0.04, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE, srt = 90)

box(bty = "l", lwd = 0.5)

straightarrow(c(-0.02, 0.32), c(0.02, 0.32), arr.pos = 1, lty = 1, lwd = 0.5, arr.length = 0.07, arr.width = 0.07, xpd = TRUE)
segments(-0.02, 0.32, 0.02, 0.32, xpd = TRUE)
text(-0.09, 0.32, "\u2640", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(-0.04, 0.32, "\u2642", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(-0.063, 0.32, "+", cex = 7/par("ps")/par("cex"), family = "Arial", xpd = TRUE)
text(0.09, 0.32, "DGRP", cex = 7/par("ps")/par("cex"), xpd = TRUE)

xy.cor <- cor.test((eff$OSfemale.Eff + eff$OSmale.Eff)/2, eff$DGRP.Eff)
cor.est <- formatC(xy.cor$estimate, format = "f", digits = 2)
cor.p <- formatC(xy.cor$p.value, format = "e", digits = 2)
print(cor.p)

text(0, 0.28, parse(text = paste("paste(italic(r), \" = \", ", cor.est, ", \" (\", italic(P), \" = \", ", cor.p,  ", \")\" )")), cex = 7/par("ps")/par("cex"), xpd = TRUE)

text(grconvertX(0.05 + file.width/25.4/4*3, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("d")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


dev.off()