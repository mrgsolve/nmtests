nmtest:
	Rscript -e ".libPaths('/data/Rlibs/'); Sys.setenv(RSTUDIO_PANDOC = '/usr/lib/rstudio-server/bin/pandoc'); rmarkdown::render('nmtest7.R')"
pdf:
	Rscript -e ".libPaths('/data/Rlibs/'); Sys.setenv(RSTUDIO_PANDOC = '/usr/lib/rstudio-server/bin/pandoc'); rmarkdown::render('nmtest7.R',output_format='pdf_document')"