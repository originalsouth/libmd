include ../../Makeheader

regressiontester: regressiontester.cc rtgroups/* ../../libmd.cc ../../libmd.h ../../libmd-src/*
	$(CC) $(CCFLAGS) -o regressiontester regressiontester.cc

clean:
	rm regressiontester
