#!/usr/bin/env bash
#SBATCH --job-name=0-5
#SBATCH --array=1-50
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4096M
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err
#SBATCH --account=rrg-jeon-ac

njob=$SLURM_ARRAY_TASK_ID
cd job-$njob

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# MUSIC
mv music* music_input_2

# set jets seed
sed -i "s/SET_SEED/$njob/g" music_input_2

# hydro evolution
./mpihydro music_input_2 1>mode_2.log 2>mode_2.err
./sweeper.sh results

if grep -Fxq "All cells frozen out. Exiting." results/mode_2.log
then

# iSS
tar -zxf EOS.tar.gz
cp results/iSS_parameters.dat parameters.dat
mv results/surface_eps_0.1928.dat results/surface.dat
./iSS.e

#Add Jet_Hadrons
mv OSCAR.DAT OSCAR_soft.dat
cp music_EOS/pdg-urqmd_v3.3+.dat .
cp results/hadrons_list.dat .
./append_jet_hadrons.e
rm -fr OSCAR_soft.dat
rm pdg-urqmd_v3.3+.dat
rm hadrons_list.dat

# UrQMD
./osc2u.e < OSCAR.DAT
mv fort.14 OSCAR.input
rm -fr OSCAR.DAT
./runqmd.sh

# Toolkit
cd hadronic_afterburner_toolkit
rm -fr results
mkdir results
mv ../particle_list.dat results/particle_list.dat



    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=211 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=0 rap_min=-0.5 rap_max=0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=211 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-0.5 rap_max=0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=211 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-1 rap_max=1 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=211 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-2.5 rap_max=2.5 >> output.log

    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=-211 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=0 rap_min=-0.5 rap_max=0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=-211 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-0.5 rap_max=0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=-211 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-1 rap_max=1 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=-211 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-2.5 rap_max=2.5 >> output.log

    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=321 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=0 rap_min=-0.5 rap_max=0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=321 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-0.5 rap_max=0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=321 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-1 rap_max=1 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=321 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-2.5 rap_max=2.5 >> output.log

    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=-321 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=0 rap_min=-0.5 rap_max=0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=-321 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-0.5 rap_max=0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=-321 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-1 rap_max=1 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=-321 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-2.5 rap_max=2.5 >> output.log

    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=2212 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=0 rap_min=-0.5 rap_max=0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=2212 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-0.5 rap_max=0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=2212 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-1 rap_max=1 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=2212 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-2.5 rap_max=2.5 >> output.log

    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=-2212 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=0 rap_min=-0.5 rap_max=0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=-2212 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-0.5 rap_max=0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=-2212 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-1 rap_max=1 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=-2212 resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-2.5 rap_max=2.5 >> output.log

    # charged hadrons
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-0.5 rap_max=0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-1.0 rap_max=1.0 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-2.5 rap_max=-0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=0.5 rap_max=2.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-2.5 rap_max=2.5 >> output.log

    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=1 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=1.5 rap_max=2.5 >> output.log
############# block for rn calculation

    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-4.5 rap_max=-4.0 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-4.0 rap_max=-3.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-3.5 rap_max=-3.0 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-3.0 rap_max=-2.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-2.5 rap_max=-2.0 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-2.0 rap_max=-1.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-1.5 rap_max=-1.0 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-1.0 rap_max=-0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-0.5 rap_max=0.0 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=0.0 rap_max=0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=0.5 rap_max=1.0 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=1.0 rap_max=1.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=1.5 rap_max=2.0 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=2.0 rap_max=2.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=2.5 rap_max=3.0 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=3.0 rap_max=3.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=3.5 rap_max=4.0 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=4.0 rap_max=4.5 >> output.log

######################

mv results/particle_list.dat ../../UrQMD_results/particle_list_$njob.dat
mv results ../../spvn_results/event_$njob
cd ..

else

echo "Cannot continue, hydro job did not finish!"

fi

