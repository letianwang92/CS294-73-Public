doxy:
	doxygen Doxyfile
	rm -f doxygenDocumentation.html
	ln -s ./html/index.html doxygenDocumentation.html

view:
	xdg-open doxygenDocumentation.html

clean:
	rm -r ./latex/ ./html/
