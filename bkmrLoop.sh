#!/bin/bash -l

for size in 300 600; do
    if [[ $size == "300" ]]; then
        time_limit="3:00:00"
        mem_gb="5"
    else
        time_limit="6:00:00"
        mem_gb="10"
    fi

    for lod_quantile in 0.25 0.5 0.75; do
        for exposure_dist in lnorm unif; do
            for h_dist in linear nonlinear; do
                sbatch --time="$time_limit" --mem="${mem_gb}gb" bkmrLoD.sh $size $lod_quantile $exposure_dist $h_dist
            done
        done
    done
done

