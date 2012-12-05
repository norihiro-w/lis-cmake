all: en ja

en: pdf-en

pdf-en: dvi-en
	dvipdfmx lis-manual-en.dvi 

dvi-en:	
	latex lis-manual-en.tex
	latex lis-manual-en.tex

ja: pdf-ja

pdf-ja: dvi-ja
	dvipdfmx lis-manual-ja.dvi 

dvi-ja:	
	platex lis-manual-ja.tex
	platex lis-manual-ja.tex

clean: 
	-$(RM) *.dvi *.aux *.toc *.log

distclean: clean
	-$(RM) *.pdf
