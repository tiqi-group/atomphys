.PHONY: format format_inplace lint test

format:
	poetry run ruff format --check atomphys

lint:
	poetry run ruff check atomphys

test:
	poetry run pytest tests/