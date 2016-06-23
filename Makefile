all:

lint:
	@echo "    Linting flooddrake codebase"
	@flake8 flooddrake
	@echo "    Linting flooddrake test suite"
	@flake8 tests
	@echo "    Linting flooddrake demo suite"
	@flake8 examples

test:
	@echo "    Running all tests"
	@py.test tests $(PYTEST_ARGS)
