#!/bin/bash -l
for size in 400; do
    if [[ $size == "400" ]]; then
        time_limit="8:00:00"
        mem_gb="5"
    else
        time_limit="36:00:00"
        mem_gb="5"
    fi
    for h_func in 1 4; do
        for lod_quantile in 0.1 0.4 0.7; do
            for exposure_dist in norm; do
                for mean_offset in 0; do
                    for scale in 0.75; do
                        for correlation in ind within across; do
                            sbatch --time="$time_limit" --mem="${mem_gb}gb" bkmrLoD.sh $size $lod_quantile $exposure_dist $mean_offset $h_func $scale $correlation
                        done
                    done
                done
            done
        done
    done
done

