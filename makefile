install:
	cd bin/utils; make install

RUNDIR = run

distclean:
	cd bin/utils; make distclean
	rm -rf run* *.diff *~


rundir:
	mkdir ${RUNDIR}
	cd ${RUNDIR}; ln -s ../bin/*.py .; mkdir results

TESTDIR = run_test

test:
	@rm -f *.diff
	-@make test_los

test_los:
	@echo "test_los_rundir..." > test_los.diff
	make test_los_rundir
	@echo "test_los_run..." >> test_los.diff
	make test_los_run
	@echo "test_los_check..." >> test_los.diff
	make test_los_check

test_los_rundir:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR}

test_los_run:
	cd ${TESTDIR}; python pyLOS_vectorized.py > runlog

test_los_check:
	gunzip -c output/result_0.txt.gz > ${TESTDIR}/results/result_0.ref
	diff ${TESTDIR}/results/result_0.txt ${TESTDIR}/results/result_0.ref \
		> test_los.diff
	ls -l test_los.diff
