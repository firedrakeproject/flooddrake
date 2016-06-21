all: modules

modules:
	@echo "    Installing modules"
	@python setup.py install

lint:
	@echo "    Linting flooddrake codebase"
	@flake8 flooddrake
	@echo "    Linting flooddrake test suite"
	@flake8 tests

test: modules
	@echo "    Running all tests"
	@py.test tests $(PYTEST_ARGS)



