nba2_az325.pdf: nba2_az325.Rnw
	R --vanilla -e "library(knitr); knit2pdf('nba2_az325.Rnw');"
	pdflatex nba2_az325.tex

clean:
	rm -f nba2_az325.aux nba2_az325.log nba2_az325.out nba2_az325.tex
	rm -rf figure
	rm -f *~
