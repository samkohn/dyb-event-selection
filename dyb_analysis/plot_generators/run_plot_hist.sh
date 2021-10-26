module load root
conda activate dyb
in_path=$1
out_path=$2
echo "Input path: " $in_path
echo "Output path: " $out_path
sleep 10
for eh in 1 2; do for ad in 1 2; do
    base_path=$in_path/EH${eh}/double_coincidences_ad${ad}
    infile=$base_path/double_coincs_EH${eh}_AD${ad}.root
    if [ ! -f $infile ]; then
        hadd -f $infile $base_path/double_coincs_EH${eh}_AD${ad}*.root
    fi
    python plot_hist.py \
        $infile double_coincidences \
        --ehad ${eh} ${ad} --logz --colz --rebin2d 5 5 \
        --xlabel "Prompt energy [MeV]" --ylabel "Delayed energy [MeV]" \
        -o $out_path/double_coincs_EH${eh}_AD${ad}.pdf
done; done

for eh in 3; do for ad in 1 2 3 4; do
    base_path=$in_path/EH${eh}/double_coincidences_ad${ad}
    infile=$base_path/double_coincs_EH${eh}_AD${ad}.root
    if [ ! -f $infile ]; then
        hadd -f $infile $base_path/double_coincs_EH${eh}_AD${ad}*.root
    fi
    python plot_hist.py \
        $infile double_coincidences \
        --ehad ${eh} ${ad} --logz --colz --rebin2d 5 5 \
        --xlabel "Prompt energy [MeV]" --ylabel "Delayed energy [MeV]" \
        -o $out_path/double_coincs_EH${eh}_AD${ad}.pdf
done; done

for eh in 1 2; do for ad in 1 2; do
    base_path=$in_path/EH${eh}/post_DT_cut_ad${ad}
    infile=$base_path/post_DT_cut_EH${eh}_AD${ad}.root
    if [ ! -f $infile ]; then
        hadd -f $infile $base_path/post_DT_cut_EH${eh}_AD${ad}*.root
    fi
    python plot_hist.py \
        $infile double_coincidences \
        --ehad ${eh} ${ad} --logz --colz --rebin2d 5 5 \
        --xlabel "Prompt energy [MeV]" --ylabel "Delayed energy [MeV]" \
        -o $out_path/post_DT_cut_EH${eh}_AD${ad}.pdf
done; done

for eh in 3; do for ad in 1 2 3 4; do
    base_path=$in_path/EH${eh}/post_DT_cut_ad${ad}
    infile=$base_path/post_DT_cut_EH${eh}_AD${ad}.root
    if [ ! -f $infile ]; then
        hadd -f $infile $base_path/post_DT_cut_EH${eh}_AD${ad}*.root
    fi
    python plot_hist.py \
        $infile double_coincidences \
        --ehad ${eh} ${ad} --logz --colz --rebin2d 5 5 \
        --xlabel "Prompt energy [MeV]" --ylabel "Delayed energy [MeV]" \
        -o $out_path/post_DT_cut_EH${eh}_AD${ad}.pdf
done; done

for eh in 1 2; do for ad in 1 2; do
    base_path=$in_path/EH${eh}/dr_vs_dt_ad${ad}
    infile=$base_path/dr_vs_dt_EH${eh}_AD${ad}.root
    if [ ! -f $infile ]; then
        hadd -f $infile $base_path/dr_vs_dt_EH${eh}_AD${ad}*.root
    fi
    python plot_hist.py \
        $infile dr_vs_dt \
        --ehad ${eh} ${ad} --logz --colz --rebin2d 2 1 \
        --xlabel 'Coincidence time [#mus]' --ylabel 'Coincidence distance [m]' \
        --xrange 0 600 --yrange 0 1 \
        -o $out_path/dr_vs_dt_EH${eh}_AD${ad}.pdf
done; done

for eh in 3; do for ad in 1 2 3 4; do
    base_path=$in_path/EH${eh}/dr_vs_dt_ad${ad}
    infile=$base_path/dr_vs_dt_EH${eh}_AD${ad}.root
    if [ ! -f $infile ]; then
        hadd -f $infile $base_path/dr_vs_dt_EH${eh}_AD${ad}*.root
    fi
    python plot_hist.py \
        $infile dr_vs_dt \
        --ehad ${eh} ${ad} --logz --colz --rebin2d 2 1 \
        --xlabel 'Coincidence time [#mus]' --ylabel 'Coincidence distance [m]' \
        --xrange 0 600 --yrange 0 1 \
        -o $out_path/dr_vs_dt_EH${eh}_AD${ad}.pdf
done; done

#for eh in 1 2 3
#do
    #python plot_hist.py $SCRATCH/dyb30/processed/EH${eh}/sub_ad1/sub_ad1.root ed_DT_sub \
        #--logz --colz -o ed_DT_sub_EH${eh}_AD1.pdf --ehad $eh 1 --rebin2d 1 10 \
        #--xlabel "DT value [mm]" --ylabel "Delayed energy [MeV]" --smallxnumbers
#done

