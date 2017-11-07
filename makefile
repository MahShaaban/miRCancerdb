all: build_db launch_app
	
build_db:
	R CMD BATCH build_script.R

launch_app:
	R -e "shiny::runApp('app.R', launch.browser = TRUE)"
