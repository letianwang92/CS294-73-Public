#!/bin/bash
doxygen Doxyfile
rm -f doxygenDocumentation.html
ln -s ./html/index.html doxygenDocumentation.html
