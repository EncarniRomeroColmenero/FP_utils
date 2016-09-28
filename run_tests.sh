#!/bin/bash

ERROR=0

python -m unittest discover tests/unittests
if [[ $? != 0 ]]
then
    ERROR=1
    echo 'Unit tests did NOT pass!' >&2
fi

exit $ERROR
