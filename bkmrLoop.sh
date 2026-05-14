#!/bin/bash -l
for size in 400; do
    if [[ $size == "400" ]]; then
        time_limit="2:00:00"
        mem_gb="3"
    else
        time_limit="8:00:00"
        mem_gb="3"
    fi

    for lod_quantile in 0.25 0.50 0.75; do
        for exposure_dist in lnorm unif; do
            for mean_offset in -1 0 1; do
                sbatch --time="$time_limit" --mem="${mem_gb}gb" bkmrLoD.sh $size $lod_quantile $exposure_dist $mean_offset
            done
        done
    done
done

