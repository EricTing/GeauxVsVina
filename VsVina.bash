#!/bin/bash

readonly ENV=/work/jaydy/working/InitialSdf

readonly SAMPLE_DIR=/work/jaydy/working/Traces

readonly ASTEX_DIR=/home/jaydy/work/dat/astex

readonly PROTEIN_DIR=/home/jaydy/work/protein-pdbqt

readonly LST=/work/jaydy/working/total_lst

readonly SDF_GEN=/home/jaydy/Workspace/script/Bashscripts/trace.bash # generate sdf file

readonly RUN_VINA=/home/jaydy/Workspace/script/Bashscripts/run_vina.bash

readonly GYRATION=/home/jaydy/Workspace/script/Perlscripts/box-gyration.pl

readonly VINA_BIN=/project/michal/apps/autodock_vina_1_1_2_linux_x86/bin/vina

readonly GEO_OFN=/work/jaydy/working/InitialSdf/vina_search_space.txt

readonly TRANS_DIR=/work/jaydy/working/Trans

readonly VINA_OUT=/work/jaydy/working/VinaResult

readonly RESTORE_ORDER=/home/jaydy/Workspace/GeauxVsVina/restore_order.py


echo_init_state() {
    # echo the initial state vector to the output file
    local complex=$1
    local sample_fn=${SAMPLE_DIR}/${complex}.sample
    local init_state=$(grep "^0.0 0.0" $sample_fn | head -n 1)
    local out_ofn=${ENV}/${complex}_init.txt
    echo $init_state > $out_ofn
}


binding_site_geo() {
    local complex=$1
    local native_sdf=${ASTEX_DIR}/${complex}/${complex}.sdf
    local gyro_out=$(perl $GYRATION $native_sdf)
    echo $complex $gyro_out
}


sdf_to_pdbqt() {
    # convert ligand sdf file to pdbqt file using babel
    local complex=$1
    local sdf=${ENV}/${complex}_init.sdf
    local trans_dir=${TRANS_DIR}/${complex}
    mkdir -p $trans_dir
    cp $sdf $trans_dir
    local id=${trans_dir}/$complex
    babel -isdf $sdf -opdbqt -h ${id}_0.pdbqt
    babel -ipdbqt ${id}_0.pdbqt -osdf ${id}_0.sdf
    babel -isdf ${id}_0.sdf -opdbqt -h ${id}_1.pdbqt
    babel -ipdbqt ${id}_1.pdbqt -osdf ${id}_1.sdf
    babel -isdf ${id}_1.sdf -opdbqt -h ${id}_2.pdbqt
    babel -ipdbqt ${id}_2.pdbqt -osdf ${id}_2.sdf
    babel -isdf ${id}_2.sdf -opdbqt -h ${id}_3.pdbqt
    babel -ipdbqt ${id}_3.pdbqt -osdf ${id}_3.sdf
    babel -isdf ${id}_3.sdf -opdbqt -h ${id}_4.pdbqt
    babel -ipdbqt ${id}_4.pdbqt -osdf ${id}_4.sdf
    babel -isdf ${id}_4.sdf -opdbqt -h ${id}_5.pdbqt
    babel -ipdbqt ${id}_5.pdbqt -osdf ${id}_5.sdf
    babel -isdf ${id}_5.sdf -opdbqt -h ${id}_6.pdbqt
}

diff_pdbqt() {
    local complex=$1
    local trans_dir=${TRANS_DIR}/${complex}
    mkdir -p $trans_dir
    local id=${trans_dir}/$complex
    diff ${id}_5.pdbqt ${id}_6.pdbqt
}

check_order_before_after_vina() {
    local complex=$1
    local output_dir=${VINA_OUT}/${complex}
    local input=${ENV}/${complex}_init.pdbqt
    local before=$(grep "ATOM" $input |awk '{print $2,$3,$4,$11,$12}')
    for vina_result in `ls ${output_dir}/*out.pdbqt`
    do
	      after=$(awk '/MODEL 1/{f=1;next}/MODEL 2/{f=0}f' $vina_result | \
	                     grep ATOM | \
	                     awk '{print $2,$3,$4,$11,$12}')

	      if [ "$before" != "$after" ]; then
	          echo $vina_result "is different"
	      fi
    done
}


highest_cms_of_sampling() {
    # to find the highest cms in running AutoVina to sample this complex
    local complex=$1
    local output=${VINA_OUT}/${complex}/vina_cms.out
    local highest_cms=$(grep "predicted conformation cms value" $output | \
	                             awk '{print $5}' | sort | tail -n 1)
    echo $complex $highest_cms
}

main() {
    cd $ENV
    for complex in `cat $LST`
    do
	      # echo_init_state $complex
	      # bash $SDF_GEN $complex
	      # sdf_to_pdbqt $complex
	      # diff_pdbqt $complex
	      # check_order_before_after_vina $complex
	      # bash $RUN_VINA $complex
	      # python $RESTORE_ORDER $complex
	      highest_cms_of_sampling $complex
	      # echo $complex "...done"
    done
}


main
