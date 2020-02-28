#!/usr/bin/env python

import sys
from os import path

class color:
    """
    define colors in the terminal
    """
    purple = '\033[95m'
    cyan = '\033[96m'
    darkcyan = '\033[36m'
    blue = '\033[94m'
    green = '\033[92m'
    yellow = '\033[93m'
    red = '\033[91m'
    bold = '\033[1m'
    underline = '\033[4m'
    end = '\033[0m'

def write_analysis_spectra_and_vn_commands(script, after_burner_type):
    #pid_particle_list = ['211', '-211', '321', '-321', '2212', '-2212',
    #                     '3122', '-3122', '3312', '-3312', '3334', '-3334',
    #                     '333']
    #charged_particle_list = ['9998', '-9998', '9999']
    pid_particle_list = ['211', '-211', '321', '-321', '2212', '-2212']
    charged_particle_list = ['9999']

    read_in_mode = 1
    if after_burner_type == "JAM":
        read_in_mode = 5
    if after_burner_type == "OSCAR":
        read_in_mode = 0

    for ipart in pid_particle_list:
        script.write(
"""
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=0 rap_min=-0.5 rap_max=0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-0.5 rap_max=0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-1 rap_max=1 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-2.5 rap_max=2.5 >> output.log
""".format(read_in_mode, ipart))
    for ipart in charged_particle_list:
        script.write(
"""
    # charged hadrons
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-0.5 rap_max=0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-1.0 rap_max=1.0 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-2.5 rap_max=-0.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=0.5 rap_max=2.5 >> output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-2.5 rap_max=2.5 >> output.log
""".format(read_in_mode, ipart))

def generate_submit_script():
    print(color.purple + "\n" + "-"*80
          + "\n>>>>> generating submission script! <<<<<\n" + "-"*80
          + color.end)
    ppn = 16
    firstA = '%A'
    secondA = '%a'
    walltime = '12:00:00'
    working_folder = path.abspath('./')
    folder_name = working_folder.split('/')[-1]
    script = open("submit_full_job.sh", "w")
    script.write(
"""#!/usr/bin/env bash
#SBATCH --job-name=%s
#SBATCH --array=1-35
#SBATCH --time=%s
#SBATCH --cpus-per-task=%d
#SBATCH --mem-per-cpu=4096M
#SBATCH --output=slurm-%s_%s.out
#SBATCH --error=slurm-%s_%s.err
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
cp results/iSS_parameters.dat .
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


""" % (folder_name, walltime, ppn, firstA, secondA, firstA, secondA))

    write_analysis_spectra_and_vn_commands(script, "UrQMD")

    script.write(
"""
mv results/particle_list.dat ../../UrQMD_results/particle_list_$njob.dat
mv results ../../spvn_results/event_$njob
cd ..

else

echo "Cannot continue, hydro job did not finish!"

fi

""")
    script.close()

generate_submit_script()
