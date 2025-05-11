#!/usr/bin/env bash 

uv run pytest --cov nanocompore --cov-branch --cov-report term-missing --cov-report html -n 4 tests/
