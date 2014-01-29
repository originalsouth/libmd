tests:
	make test-orbit
	make test-curved-orbit
	make test-randomspring

test-orbit:
	make -C ./orbit

test-curved-orbit:
	make -C ./curved-orbit

test-randomspring:
	make -C ./randomspring

clean:
	make cleanorbit
	make cleancurved-orbit
	make cleanrandomspring

cleanorbit:
	make -C ./orbit clean

cleancurved-orbit:
	make -C ./curved-orbit clean

cleanrandomspring:
	make -C ./randomspring clean
