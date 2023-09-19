.PHONY: format format_inplace lint test

format:
	poetry run autopep8 --diff --recursive --exclude _*.py --exit-code atomphys/

format_inplace:
	poetry run autopep8 --in-place --recursive --exclude _*.py --exit-code atomphys/

lint:
	poetry run flake8 --config setup.cfg atomphys/

test:
	poetry run pytest tests/