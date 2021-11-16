#! /bin/bash

# Script base directory.
basedir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Generate the binaries
PYTHON=/opt/python/cp37-cp37m/bin/python3
SCRIPT=$(cat <<-END
echo "installing build dependencies"
cd /work
yum install -y libffi-devel
$PYTHON -m pip install --target=dist/site-packages -U cffi pcpp
export PYTHONPATH=/work/dist/site-packages

echo "building libalouette with \$(gcc --version | head -n 1)"
make PYTHON=$PYTHON

echo "building Python wheel"
make package PYTHON=$PYTHON
$PYTHON setup.py bdist_wheel --py-limited-api=cp37
auditwheel show dist/*.whl
auditwheel repair dist/*.whl

echo "cleaning build artifacts"
make dist-clean
mv wheelhouse/*.whl dist
rm -rf wheelhouse
chown --recursive $(id -u):$(id -g) alouette build dist lib
END
)

make clean
make dist-clean
mkdir -p dist/site-packages
docker run --mount type=bind,source=$(pwd),target=/work                        \
                   quay.io/pypa/manylinux1_x86_64 /bin/bash -c "${SCRIPT}"
