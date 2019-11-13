### Before Starting

Each new user who runs the pipelines must also first update/install some libraries in python and R.

    #On the cluster, load the python module to get access to python3 and a wide variety of libraries
    module load python3.7-anaconda/2019.07

    #Update snakemake, installing the update to user’s local directory (pipeline tested w/ v5.4.5)
    pip install snakemake -U --user
    #Update pandas, installing to user’s local directory (pipeline tested w/ v0.24)
    pip install pandas -U --user


    #Install optparse R library to a local_R directory
    mkdir ~/local_R
    module load R/3.6.1 ; export R_LIBS_USER=$HOME/local_R
    R
    #Within R console:
    install.packages("optparse")
    quit()

### Quick-start example

The config files for this pipeline are manually generated. An example config resides in the example folder.
Modify the last three entries to change the locations where logs, results, and temp files are placed.

    #Start screen or tmux to create a persistent session over ssh (I prefer screen)
    screen
    #Navigate to the location of the snakefile
    cd /nfs/turbo/path-rjhryan-turbo/lab-members/Travis/projects/snakemake_CRISPRi
    #Modify the config as described above, then
    #Perform a dry-run (-n flag for dry-run, -p flag to print shell commands)
    snakemake -n -p --snakefile Snakefile_CRISPRi --configfile example/CRISPRi_2334_example_config.json

    #If that looks OK, then run the pipeline on the cluster
    snakemake --snakefile Snakefile_CRISPRi --configfile example/CRISPRi_2334_example_config.json \
    -j 144 --cluster-config cluster_config.json \
    --cluster 'sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time} --output=%x-%j.out'


### More info

There also exist several old run configs in the run_configs subfolder, which were used on earlier versions of the pipeline. They should still work, though I haven't explicitly tested them.

Aside from the bowtie step, this pipeline uses very little in the way of resources, and could probably even be run on the login node, e.g.

    snakemake --snakefile Snakefile_CRISPRi --configfile example/CRISPRi_2334_example_config.json
