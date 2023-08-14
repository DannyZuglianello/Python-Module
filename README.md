# pylibmaterials

Python Interface to laurent/libmaterials.

To compile:

mkdir -p gitignore/ && g++ -fpic -c pyplugin.cpp -o ./gitignore/pyplugin.o && g++ -shared ./gitignore/pyplugin.o -o ./gitignore/pylibmaterials.so
