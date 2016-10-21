all:

docs:
	@echo "    generating flooddrake documentation"
	@pdflatex -interaction batchmode flooddrake_doc
	@bibtex flooddrake_doc
	@pdflatex -interaction batchmode flooddrake_doc
	@pdflatex -interaction batchmode flooddrake_doc
