all: build_db launch_app
	
build_db:
	Rscript R/build_script.R

clean:
	rm -rf tmp

launch_app:
	R -e "shiny::runApp('app.R', launch.browser = TRUE)"
