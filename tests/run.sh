#! /bin/bash
PYTHON=${PYTHON:-python3}
TESTS=${TESTS:-$(ls tests/test_*.py)}

coverage="$PYTHON -m coverage"

$coverage erase
for filename in $TESTS; do
    $coverage run --source=alouette -a -m pytest $filename
done
$coverage html
