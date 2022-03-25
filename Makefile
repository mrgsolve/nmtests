all: 
	make nmtest
	make pdf
nmtest:
	Rscript -e "Sys.setenv(RSTUDIO_PANDOC = '/usr/lib/rstudio-server/bin/pandoc'); rmarkdown::render('nmtest9.R')"
pdf:
	Rscript -e "Sys.setenv(RSTUDIO_PANDOC = '/usr/lib/rstudio-server/bin/pandoc'); rmarkdown::render('nmtest9.R',output_format='pdf_document')"