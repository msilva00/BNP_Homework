## Terminal Instructions for DPpackage installation for Mac Users:
### Check if you have fortran:
	which gfortran
	# /usr/local/bin/gfortran  ## <-- My output: I already had it but it was old

### Check the version
	gfortran -v  ## check the version
	 # it outputs a bunch of stuff, says version 4.8 … I don’t want that one
### remove the old/unwanted version	 
	rm /usr/local/bin/gfortran # get rid of this version
	brew reinstall gcc # this has newest version of gfortran and C/C++ compiler
### If you've never installed a C/C++ compiler.
	#brew install gcc ## uncomment and run if you have never installed gcc
	# mkdir ~/.R ##<— do if you don’t have directory mkdir ~/.R
### Else do the following (change Makevar from Terminal)
Be careful with the following as it may change existing Makevars/settings.

	cat << EOF >> ~/.R/Makevars
	FLIBS=-L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin16/6.3.0 -L/usr/local/gfortran/lib -lgfortran -lquadmath -lm
	EOF


## R instructions (run in R studio):
	url = "https://cran.r-project.org/src/contrib/Archive/DPpackage/DPpackage_1.1-7.4.tar.gz" # most recent archive
	pkgFile = "DPpackage_1.1-7.4.tar.gz"
	download.file(url = url, destfile = pkgFile)
	install.packages(c("ada", "ipred", "evd"))
	install.packages(pkgs=pkgFile, type="source", repo=NULL)
