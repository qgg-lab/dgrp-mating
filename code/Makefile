# ===============================================
# = make file for generating tables and figures =
# ===============================================

all: ../figures/figureRawData.pdf ../figures/figureCorr.pdf ../figures/figureSNPeff.pdf
	
FONT="Myriad Pro"

# figure for raw data distribution
# ============================================================

../figures/figureRawData.pdf: figureRawData.R ../data/OSDGRPmating.csv
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figureRawData.R ../figures/figureRawData.pdf $(FONT) ../data/OSDGRPmating.csv

../figures/figureCorr.pdf: figureCorr.R ../data/OSDGRPmating.csv ../data/cross.dgrp.csv
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figureCorr.R ../figures/figureCorr.pdf $(FONT) ../data/OSDGRPmating.csv ../data/cross.dgrp.csv

../figures/figureSNPeff.pdf: figureSNPeff.R ../data/DGRPxOSeff.csv
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figureSNPeff.R ../figures/figureSNPeff.pdf $(FONT) ../data/DGRPxOSeff.csv
