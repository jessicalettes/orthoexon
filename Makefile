help:
	@echo "lint - check code style with flake8"
	@echo "test - run tests quickly"
	@echo "coverage - check code coverage quickly"

test:
	py.test orthoexon

coverage:
	coverage run --source orthoexon --omit="*/test*" --module py.test
	coverage report --show-missing

lint:
	flake8 orthoexon
