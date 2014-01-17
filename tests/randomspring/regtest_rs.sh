#!/bin/bash

make STD=c++0x randomspring
./randomspring >/dev/null
cat sim*.pts sim*.bds | md5sum > __foo
cat __foo hash_check
diff __foo hash_check
rm __foo

