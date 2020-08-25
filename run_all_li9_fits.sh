python li9_fit.py li9_fit_low.pdf $SCRATCH/dyb30/parameters.db $SCRATCH/dyb30/processed/EH1/li9_hadded_ad?/li9_*.root --energy low --site 1 --update-db &
python li9_fit.py li9_fit_low_ntag.pdf $SCRATCH/dyb30/parameters.db $SCRATCH/dyb30/processed/EH1/li9_hadded_ad?/li9_*.root --energy low --ntag --site 1 --update-db &
python li9_fit.py li9_fit_mid.pdf $SCRATCH/dyb30/parameters.db $SCRATCH/dyb30/processed/EH1/li9_hadded_ad?/li9_*.root --energy mid --site 1 --update-db &
python li9_fit.py li9_fit_mid_ntag.pdf $SCRATCH/dyb30/parameters.db $SCRATCH/dyb30/processed/EH1/li9_hadded_ad?/li9_*.root --energy mid --ntag --site 1 --update-db &
python li9_fit.py li9_fit_high.pdf $SCRATCH/dyb30/parameters.db $SCRATCH/dyb30/processed/EH1/li9_hadded_ad?/li9_*.root --energy high --site 1 --update-db &
python li9_fit.py li9_fit_high_ntag.pdf $SCRATCH/dyb30/parameters.db $SCRATCH/dyb30/processed/EH1/li9_hadded_ad?/li9_*.root --energy high --ntag --site 1 --update-db &
wait
