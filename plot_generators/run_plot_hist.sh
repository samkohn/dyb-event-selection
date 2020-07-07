#for eh in 1 2; do for ad in 1 2; do python plot_hist.py \
    #double_coincs_EH${eh}_AD${ad}.root double_coincidences --logz --colz \
    #--rebin2d 5 5; done; done

#for ad in 1 2 3 4; do python plot_hist.py double_coincs_EH3_AD${ad}.root \
    #double_coincidences --logz --colz --rebin2d 5 5; done

for eh in 1 2; do for ad in 1 2; do python plot_hist.py \
    post_DT_cut_EH${eh}_AD${ad}.root double_coincidences --logz --colz \
    --rebin2d 5 5; done; done

for ad in 1 2 3 4; do python plot_hist.py post_DT_cut_EH3_AD${ad}.root \
    double_coincidences --logz --colz --rebin2d 5 5; done

#for eh in 1 2; do for ad in 1 2; do python plot_hist.py \
    #dr_vs_dt_EH${eh}_AD${ad}.root dr_vs_dt --logz --colz --rebin2d 2 1 \
    #--xrange 0 600 --yrange 0 1.5; done; done

#for ad in 1 2 3 4; do python plot_hist.py dr_vs_dt_EH3_AD${ad}.root dr_vs_dt \
    #--logz --colz --rebin2d 2 1 --xrange 0 600 --yrange 0 1.5; done

#for eh in 1 2 3
#do
    #python plot_hist.py $SCRATCH/dyb30/processed/EH${eh}/sub_ad1/sub_ad1.root ed_DT_sub \
        #--logz --colz -o ed_DT_sub_EH${eh}_AD1.pdf --ehad $eh 1 --rebin2d 1 10 \
        #--xlabel "DT value [mm]" --ylabel "Delayed energy [MeV]" --smallxnumbers
#done

