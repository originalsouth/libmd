tests:
	make test-regressiontester
	make test-boxshear
	make test-forcedipole

test-regressiontester:
	make -C ./regressiontester

test-boxshear:
	make -C ./boxshear

test-forcedipole:
	make -C ./forcedipole

clean:
	make cleanregressiontester
	make cleanboxshear
	make cleanforcedipole 

cleanregressiontester:
	make -C ./regressiontester clean

cleanboxshear:
	make -C ./boxshear clean

cleanforcedipole:
	make -C ./forcedipole clean

