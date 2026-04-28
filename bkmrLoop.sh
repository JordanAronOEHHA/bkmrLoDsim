#!/bin/bash -l
for size in 800; do
    if [[ $size == "400" ]]; then
        time_limit="3:00:00"
        mem_gb="3"
    else
        time_limit="16:00:00"
        mem_gb="3"
    fi

    for lod_quantile in 0.2 0.4 0.6 0.8; do
        for exposure_dist in gamma lnorm unif; do
            for h_dist in linear nonlinear; do
                sbatch --time="$time_limit" --mem="${mem_gb}gb" bkmrLoD.sh $size $lod_quantile $exposure_dist $h_dist
            done
        done
    done
done

