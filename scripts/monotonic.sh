#!/bin/bash

foo=0
while read p; do
  foo=$((foo+1));
  if [ $foo -ne $p ]
  then
  	print "BAD"
	exit 0
  fi
done < $1
