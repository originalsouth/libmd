tests:
	make test-regressiontester
	make test-boxshear
	make test-forcedipole
	make test-autodiff

test-regressiontester:
	make -C ./regressiontester

test-boxshear:
	make -C ./boxshear

test-forcedipole:
	make -C ./forcedipole

test-autodiff:
	make -C ./autodiff

clean:
	make cleanregressiontester
	make cleanboxshear
	make cleanforcedipole 
	make cleanautodiff

cleanregressiontester:
	make -C ./regressiontester clean

cleanboxshear:
	make -C ./boxshear clean

cleanforcedipole:
	make -C ./forcedipole clean

cleanautodiff:
	make -C ./autodiff clean
