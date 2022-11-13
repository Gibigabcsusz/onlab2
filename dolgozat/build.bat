ECHO OFF
set cim="onlab2"
echo "Forditas pdflatex 1"
echo "Forditas pdflatex 1\n\n\n" > output.txt
pdflatex -interaction=nonstopmode -halt-on-error %cim%.tex >> output.txt
echo "Forditas bibtex"
echo "\n\n\nForditas bibtex\n\n\n" >> output.txt
bibtex %cim% >> output.txt
echo "Forditas pdflatex 2"
echo "\n\n\nForditas pdflatex 2\n\n\n" >> output.txt
pdflatex -interaction=nonstopmode -halt-on-error %cim%.tex >> output.txt
echo "Forditas pdflatex 3"
echo "\n\n\nForditas pdflatex 3\n\n\n" >> output.txt
pdflatex -interaction=nonstopmode -halt-on-error %cim%.tex >> output.txt

if exist %cim%.aux rm %cim%.aux
if exist %cim%.bbl rm %cim%.bbl
if exist %cim%.blg rm %cim%.blg
if exist %cim%.out rm %cim%.out
if exist %cim%.log rm %cim%.log
if exist %cim%.toc rm %cim%.toc
if exist %cim%.xmpdata rm %cim%.xmpdata
if exist pdfa.xmpi rm pdfa.xmpi
if exist tartalom\*.aux rm tartalom\*.aux

rem rm %cim%.aux %cim%.bbl %cim%.blg %cim%.out %cim%.log %cim%.toc
rem rm tartalom/*.aux
