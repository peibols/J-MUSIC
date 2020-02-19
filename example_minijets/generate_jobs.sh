#!/bin/bash

music_folder=/home/peibols/projects/rrg-jeon-ac/peibols/source/music_stuff
init_folder=/home/peibols/projects/rrg-jeon-ac/group_writable/IP-Glasma_Updated_2D_events/2.76TeV/0-5

iss_folder=/home/peibols/projects/rrg-jeon-ac/peibols/iSS-master
urqmd_folder=/home/peibols/projects/rrg-jeon-ac/peibols/urqmd
toolkit_folder=/home/peibols/projects/rrg-jeon-ac/peibols/source/urqmd_analysis/toolkit
append_hadron_folder=/home/mayanks/projects/rrg-jeon-ac/mayanks/add_jet_hadrons

mkdir UrQMD_results
mkdir spvn_results

for ii in {1..150}
  do
    mkdir job-$ii

    # MUSIC
    cp $music_folder/mpihydro job-$ii
    cp $music_folder/sweeper.sh job-$ii
    cp $music_folder/eps_freeze_list_hotqcd.dat job-$ii
    cp $music_folder/setup_pythia.cmnd job-$ii
    cp -r $music_folder/EOS job-$ii
    cp $music_folder/music_input_minijets job-$ii
    cp $init_folder/binary_collision_locations_$ii.dat job-$ii/binary_collision_locations.dat
    cp $init_folder/u_field_$ii.dat job-$ii/u_field.dat
    
    # iSS
    cp $iss_folder/iSS.e job-$ii/
    cp $iss_folder/EOS.tar.gz job-$ii/
    cp -r $iss_folder/iSS_tables job-$ii/
    cp $iss_folder/iSS_parameters.dat job-$ii/
    cp $urqmd_folder/OSCAR_header.txt job-$ii/
    cp $append_hadron_folder/append_jet_hadrons.e job-$ii/

    # UrQMD
    cp $urqmd_folder/uqmd.burner job-$ii/
    cp $urqmd_folder/urqmd.e job-$ii/
    cp $urqmd_folder/runqmd.sh job-$ii/
    cp $urqmd_folder/osc2u.e job-$ii/
      
    # Toolkit
    cp -r $toolkit_folder/hadronic_afterburner_toolkit job-$ii/

    echo "Created folder job-$ii"
  done

  ./generate_full_inputfile.py
