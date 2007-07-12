#!/bin/bash

hc_findplate -B data -D enes -O plates.id

hc_ptrot -B data -D enes -P plates.id -O unitrots_new.coeffs 

hcplates -L point.j -T plates_ids.ixz -U unitrots_new.coeffs -Op mypoles.out -Ov myvels.out -P parameter_file.default

pmyvels_simple


