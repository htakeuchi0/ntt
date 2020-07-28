# commands
TEX = uplatex
MENDEX = upmendex
MENDEXFLAGS = -s
MENDEXFLAGS_JAPANESE = -g
DVIPDF = dvipdfmx

# pdfview
PDFVIEW = evince

# source
SRC = ntt
SRC_TEX = $(SRC).tex
SRC_DVI = $(SRC).dvi
SRC_PDF = $(SRC).pdf
SRC_IST = $(SRC).ist

# extensions
EXTS := .aux .log .toc .dvi .pdf .lof .lot .ilg .ind .idx

.PHONY: all clean

all: show

$(SRC_DVI): $(SRC_TEX)
	$(TEX) $(SRC_TEX)
	$(MENDEX) $(MENDEXFLAGS) $(SRC_IST) $(MENDEXFLAGS_JAPANESE) $(SRC)
	$(TEX) $(SRC_TEX)

$(SRC_PDF): $(SRC_DVI)
	$(DVIPDF) $(SRC_DVI)

show: $(SRC_PDF)
	$(PDFVIEW) $(SRC_PDF) &

clean:
	@for ext in $(EXTS); do \
		rm -f $(SRC)$$ext; \
	done

