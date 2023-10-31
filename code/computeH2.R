# ================
# = calculate H2 =
# ================

args <- commandArgs(TRUE) # args <- c("../data/OSDGRPmating.csv")
library("RColorBrewer")
library("lmerTest")
library("optimx")

# read raw data
# ============================================================
os.data <- read.csv(args[1], header = TRUE, as.is = TRUE, na.strings = c("NA", ""))
os.data <- subset(os.data, !is.na(series) & n.female == 5 & n.male == 10)

# make calculations
os.data <- within(os.data, {female.line <- factor(female.line); male.line <- factor(male.line); ratio.30m <- as.numeric(ratio.30m); ratio.30m2 <- ifelse(is.na(n.mate.30m), 0, n.mate.30m)/5; n.mate.1hr <- as.numeric(n.mate.1hr); ratio.1hr2 <- ifelse(is.na(n.mate.1hr), 0, n.mate.1hr)/5 })

# replace ratio.30m and ratio.1hr in the original with the new calculations
# there were some errors
os.data <- within(os.data, {ratio.30m <- ratio.30m2; ratio.1hr <- ratio.1hr2;})
dgrp.self <- subset(os.data, as.character(female.line) == as.character(male.line) & as.character(female.line) != "OS")
os.female <- subset(os.data, as.character(female.line) == "OS" & as.character(male.line) != "OS")
os.male <- subset(os.data, as.character(male.line) == "OS" & as.character(female.line) != "OS")

# fit mixed models
# ============================================================

mixed.res <- NULL
os.male.fit <- lmer(ratio.30m ~ (1|female.line), data = os.male)
mixed.res <- rbind(mixed.res,
  c("dgrp x os", as.data.frame(VarCorr(os.male.fit))[, 4],
  as.data.frame(VarCorr(os.male.fit))[1, 4]/sum(as.data.frame(VarCorr(os.male.fit))[, 4]),
  rand(os.male.fit)$Pr[2]))

os.female.fit <- lmer(ratio.1hr ~ (1|male.line), data = os.female)
mixed.res <- rbind(mixed.res,
  c("os x dgrp", as.data.frame(VarCorr(os.female.fit))[, 4],
  as.data.frame(VarCorr(os.female.fit))[1, 4]/sum(as.data.frame(VarCorr(os.female.fit))[, 4]),
  rand(os.female.fit)$Pr[2]))

dgrp.self.fit <- lmer(ratio.1hr ~ (1|male.line), data = dgrp.self)
mixed.res <- rbind(mixed.res,
  c("dgrp x dgrp", as.data.frame(VarCorr(dgrp.self.fit))[, 4],
  as.data.frame(VarCorr(dgrp.self.fit))[1, 4]/sum(as.data.frame(VarCorr(dgrp.self.fit))[, 4]),
  rand(dgrp.self.fit)$Pr[2]))

cat(mixed.res)

dev.off()