#!/bin/bash

#Run within a directory containing setup.py

python setup.py sdist bdist_wheel && python /Users/jonathanrubin/Library/Python/2.7/lib/python/site-packages/twine upload dist/*
