.PHONY: format format_inplace lint test

format:
	uv run ruff format --check atomphys

lint:
	uv run ruff check atomphys

test:
	uv run pytest tests/
