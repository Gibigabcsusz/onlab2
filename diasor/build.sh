#!/bin/sh
cim="diasor"
bib=1
echo "Fordítás pdflatex 1"
echo "Fordítás pdflatex 1\n\n\n" > output.txt
pdflatex -interaction=nonstopmode -halt-on-error $cim.tex >> output.txt
if [ $? = 0 ]; then
    if [ $bib != 0 ]; then
        echo "Fordítás bibtex"
        echo "\n\n\nFordítás bibtex\n\n\n" >> output.txt
        bibtex $cim >> output.txt
        echo "Fordítás pdflatex 2"
        echo "\n\n\nFordítás pdflatex 2\n\n\n" >> output.txt
        pdflatex -interaction=nonstopmode -halt-on-error $cim.tex >> output.txt
    fi
    echo "Fordítás pdflatex 3"
    echo "\n\n\nFordítás pdflatex 3\n\n\n" >> output.txt
    pdflatex -interaction=nonstopmode -halt-on-error $cim.tex >> output.txt
else
    echo "Error:\n\n"
    cat output.txt
fi


if [ -e $cim.nav ]; then
	rm $cim.nav
fi
if [ -e $cim.aux ]; then
	rm $cim.aux
fi
if [ -e $cim.bbl ]; then
	rm $cim.bbl
fi
if [ -e $cim.blg ]; then
	rm $cim.blg
fi
if [ -e $cim.out ]; then
	rm $cim.out
fi
if [ -e $cim.log ]; then
	rm $cim.log
fi
if [ -e $cim.toc ]; then
	rm $cim.toc
fi
if [ -e $cim.xmpdata ]; then
	rm $cim.xmpdata
fi
if [ -e pdfa.xmpi ]; then
	rm pdfa.xmpi
fi
if [ -e tartalom/*.aux ]; then
	rm tartalom/*.aux
fi

#rm $cim.aux $cim.bbl $cim.blg $cim.out $cim.log $cim.toc
#rm tartalom/*.aux
