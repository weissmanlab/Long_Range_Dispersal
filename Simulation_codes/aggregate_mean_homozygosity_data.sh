#!/bin/bash

echo -e "c alpha coal mu x0 hom homLow homHigh" > aggregated_mh_data.txt
cat  ./alpha_value*/dist*/mean_homozygosity* >> aggregated_mh_data.txt

