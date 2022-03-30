all: 
	make dose

.PHONY: dosing
dosing:
	Rscript -e "rmarkdown::render('dosing/dosing-vignette.R')"
check-dosing:
	Rscript "dosing/dosing-check.R"
tv: 
	Rscript -e "rmarkdown::render('time-varying/time-varying-vignette.R')"
check-tv: 
	Rscript "time-varying/time-varying-check.R"
pdf:
	Rscript -e "Sys.setenv(RSTUDIO_PANDOC = '/usr/lib/rstudio-server/bin/pandoc'); rmarkdown::render('nmtest9.R',output_format='pdf_document')"
confirm:
	Rscript nmtest9-check.R
