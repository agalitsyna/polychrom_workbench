lifetimes_choices = [100] #[50, 100, 200] # [100]
# separations_choices = [20, 50, 100, 250, 300, 500, 750, 1000, 1250, 2500, 20000]
### log-space with 3 + log-space with 1.5:
separations_choices = [3, 11, 33, 100, 300, 900] + [225, 150, 66,  44]
n_targeted_lefs_choices = [20, 50, 100]  + [33, 300] + [337, 225, 150, 66,  44,  29,  19,  13,   8]

sea_rises = expand(
    "data/traj1d/sea-rises.L_{lifetimes}.S_{separations}.T_{n_targeted_lefs}.notads_0_.h5py",
    lifetimes=lifetimes_choices,
    separations=separations_choices,
    n_targeted_lefs=n_targeted_lefs_choices,
)

lifetimes_choices = [100] #[50, 100, 200] # [100]
separations_choices = [3, 11, 33, 100, 300, 900] + [225, 150, 66,  44]  #[50, 100, 200, 500] # [100]
enrichments_choices = [1.2, 2, 5, 10, 20, 50] #[0, 2, 3, 5, 10, 15, 20, 25, 50, 75, 100, 200, 250, 500, 750, 1000]

fountain_takeover = expand(
    "data/traj1d/fountain-takeover.L_{lifetimes}.ST_{separations}.E_{enrichments}.notads_0_.h5py",
    lifetimes=lifetimes_choices,
    separations=separations_choices,
    enrichments=enrichments_choices,
)

lifetimes_choices = [100] #[50, 100, 200] # [100]
separations_choices = [3, 11, 33, 100, 300, 900] + [225, 150, 66,  44]  #[50, 100, 200, 500] # [100]
enrichments_choices = [1.2, 2, 5, 10, 20, 50] #[0, 2, 3, 5, 10, 15, 20, 25, 50, 75, 100, 200, 250, 500, 750, 1000]
ctcf_probs_choices = [0.01, 0.05, 0.1, 0.2, 0.25, 0.8] #, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 1.0]

growing_boundaries = expand(
    "data/traj1d/growing-boundaries.L_{lifetimes}.S_{separations}.E_{enrichments}.tads_{ctcf_probs}_.h5py",
    lifetimes=lifetimes_choices,
    separations=separations_choices,
    enrichments=enrichments_choices,
    ctcf_probs=ctcf_probs_choices,
)

rule all:
    input:
        lambda wildcards: growing_boundaries + sea_rises + fountain_takeover

rule simulate_1d:
    output:
        file="data/traj1d/{mode}.L_{lifetime}.{separation_mode}_{separation}.{targeting_mode}_{targeting_param}.{tad_mode}_{tad_param}_.h5py",
    run:
        if "sea" in wildcards.mode:
            shell("""
                python 00_1D_targeted_extrusion.py --output {output.file} \
                        -L {wildcards.lifetime} -S {wildcards.separation} -T {wildcards.targeting_param} \
                        --loading-sites 250,750,1250,1750,2250 -N 2500 -R 20 -D 101000 --step 10
                """)
        elif "fountain" in wildcards.mode:
            shell("""
                python 00_1D_targeted_extrusion.py --output {output.file} \
                    -L {wildcards.lifetime} --separation-total {wildcards.separation} -E {wildcards.targeting_param} \
                    --loading-sites 250,750,1250,1750,2250 -N 2500 -R 20 -D 101000 --step 10
                """)
        elif "boundaries" in wildcards.mode:
            shell("""
		python 00_1D_targeted_extrusion.py --output {output.file} \
                    -L {wildcards.lifetime} -S {wildcards.separation} -E {wildcards.targeting_param} \
                    --ctcf-capture-probability {wildcards.tad_param} --ctcf-release-probability 0.005 \
                    --tads 100,400,600,900,1100,1400,1600,1900,2100,2400 \
                    --loading-sites 250,750,1250,1750,2250 -N 2500 -R 20 -D 101000 --step 10
                """)
        else:
            print("""Unlisted mode requested... Nothing to be done """)

# rule simulate_3d:
#     input:
#         file="data/traj1d/{mode}.L_{lifetime}.{separation_mode}_{separation}.{targeting_mode}_{targeting_param}.{tad_mode}_{tad_param}_.h5py",
#     output:
#         file="data/traj3d/{mode}.L_{lifetime}.{separation_mode}_{separation}.{targeting_mode}_{targeting_param}.{tad_mode}_{tad_param}_/blocks_9900-9999.h5",
#     shell:
#         "python 01_3D_simulations.py -v -f -G 0 -i {input.file} -o {output.file}"
#
# rule extract_maps:
#     input:
#         file="data/traj3d/{mode}.L_{lifetime}.{separation_mode}_{separation}.{targeting_mode}_{targeting_param}.{tad_mode}_{tad_param}_/blocks_9900-9999.h5",
#     output:
#         file="data/maps2d/{mode}.L_{lifetime}.{separation_mode}_{separation}.{targeting_mode}_{targeting_param}.{tad_mode}_{tad_param}_.npy",
#     shell:
#         "python 02_2D_contactmap.py `basename {input.file}` -o {output.file}"

# for file in data/traj3d/*; do python 02_2D_contactmap.py $file ${file/traj3d/maps2d}.npy; done 