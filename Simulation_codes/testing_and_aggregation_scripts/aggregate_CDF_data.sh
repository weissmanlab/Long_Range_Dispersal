#!/bin/bash

echo -e "c alpha coal T x0 cdf cdfLow cdfHigh" > aggregated_cdf_data.txt
cat  ./alpha_value*/dist*/CDF_of_coalescence_times* >> aggregated_cdf_data.txt
