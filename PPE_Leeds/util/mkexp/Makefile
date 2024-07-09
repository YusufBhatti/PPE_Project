prefix = /usr/local

all:

check:

clean:
	$(RM) -r build dist MANIFEST *.pyc

distclean: clean

doc:

dist: __FORCE__
	python setup.py sdist

install: all
	python setup.py install --prefix=$(prefix)

__FORCE__:
