tests:
	make test-regressiontester
	make test-boxshear
	make test-forcedipole
	make test-autodiff
	make test-geodesics

test-regressiontester:
	make -C ./regressiontester

test-boxshear:
	make -C ./boxshear

test-forcedipole:
	make -C ./forcedipole

test-autodiff:
	make -C ./autodiff

test-geodesics:
	make -C ./geodesics

clean:
	make cleanregressiontester
	make cleanboxshear
	make cleanforcedipole 
	make cleanautodiff
	make cleangeodesics

cleanregressiontester:
	make -C ./regressiontester clean

cleanboxshear:
	make -C ./boxshear clean

cleanforcedipole:
	make -C ./forcedipole clean

cleanautodiff:
	make -C ./autodiff clean 

cleangeodesics:
	make -C ./geodesics clean
