#!/bin/bash

#Command line syntax: sh gitPullPush.sh "commit comment here"

git pull
git add .
git commit -m "$@"
git push
