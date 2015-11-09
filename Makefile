help:
	@echo "lint - check code style with flake8"
	@echo "test - run tests quickly"
	@echo "coverage - check code coverage quickly"

test:
	py.test

coverage:
	coverage run --source orthoexon --omit=tests --module py.test

lint:
	flake8 orthoexon
