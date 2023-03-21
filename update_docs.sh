#!/bin/bash
rm docs/source/modules.rst
sphinx-apidoc ./ -o docs/source

shopt -s extglob

cd docs
rm -rf -- !("build"|"source")
cd ../

make clean
make html

shopt -s dotglob
mv docs/build/html/* docs/