run101:
	Rscript -e ".libPaths('~/Rlibs/'); Sys.setenv(RSTUDIO_PANDOC = '/usr/lib/rstudio-server/bin/pandoc'); rmarkdown::render('nmtest.R')"
	
run2:
	Rscript -e ".libPaths('~/Rlibs/'); Sys.setenv(RSTUDIO_PANDOC = '/usr/lib/rstudio-server/bin/pandoc'); rmarkdown::render('nmtest2.R')"
	