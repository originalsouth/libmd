tests:
	make test-orbit
	make test-curved-orbit
	make test-randomspring
	make test-regressiontester
	make test-boxshear
	make test-springpluscoulomb

test-orbit:
	make -C ./orbit

test-curved-orbit:
	make -C ./curved-orbit

test-randomspring:
	make -C ./randomspring

test-regressiontester:
	make -C ./regressiontester

test-boxshear:
	make -C ./boxshear

test-springpluscoulomb:
	make -C ./springpluscoulomb

clean:
	make cleanorbit
	make cleancurved-orbit
	make cleanrandomspring
	make cleanregressiontester
	make cleanboxshear
	make cleanspringpluscoulomb

cleanorbit:
	make -C ./orbit clean

cleancurved-orbit:
	make -C ./curved-orbit clean

cleanrandomspring:
	make -C ./randomspring clean

cleanregressiontester:
	make -C ./regressiontester clean

cleanboxshear:
	make -C ./boxshear clean

cleanspringpluscoulomb:
	make -C ./springpluscoulomb clean
