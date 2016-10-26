all:

docs:
	@echo "    generating Flooddrake documentation"
	@pdflatex -interaction batchmode flooddrake_doc
	@bibtex flooddrake_doc
	@pdflatex -interaction batchmode flooddrake_doc
	@pdflatex -interaction batchmode flooddrake_doc
