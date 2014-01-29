tests:
	make test-orbit
	make test-curved-orbit
	make test-randomspring
	make test-regressiontester

test-orbit:
	make -C ./orbit

test-curved-orbit:
	make -C ./curved-orbit

test-randomspring:
	make -C ./randomspring

test-regressiontester:
	make -C ./regressiontester

clean:
	make cleanorbit
	make cleancurved-orbit
	make cleanrandomspring
	make cleanregressiontester

cleanorbit:
	make -C ./orbit clean

cleancurved-orbit:
	make -C ./curved-orbit clean

cleanrandomspring:
	make -C ./randomspring clean

cleanregressiontester:
	make -C ./regressiontester clean
