install: createRay.so

createRay.so: createRay.pyx
	python setup.py build_ext --inplace

distclean:
	rm -f createRay.so createRay.c *.pyc
