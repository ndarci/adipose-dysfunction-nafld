#!/bin/bash


direc=$1

if (("$direc" == "f")); then
	echo foo
else
	echo error
fi
