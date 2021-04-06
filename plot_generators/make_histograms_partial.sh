module load root
conda activate dyb
in_path=$1
out_path=$2
digit=$3
echo "Input path: " $in_path
echo "Output path: " $out_path
echo "digit: " $digit
sleep 10
python plot_double_coincidences.py \
"${in_path}/EH1/hadded_ad1/out_ad1_${digit}????.root" \
${out_path}/EH1/double_coincidences_ad1/double_coincs_EH1_AD1_${digit}0000s --root &
sleep 1
python plot_double_coincidences.py \
"${in_path}/EH1/hadded_ad2/out_ad2_${digit}????.root" \
${out_path}/EH1/double_coincidences_ad2/double_coincs_EH1_AD2_${digit}0000s --root &
sleep 1
python plot_double_coincidences.py \
"${in_path}/EH2/hadded_ad1/out_ad1_${digit}????.root" \
${out_path}/EH2/double_coincidences_ad1/double_coincs_EH2_AD1_${digit}0000s --root &
sleep 1
python plot_double_coincidences.py \
"${in_path}/EH2/hadded_ad2/out_ad2_${digit}????.root" \
${out_path}/EH2/double_coincidences_ad2/double_coincs_EH2_AD2_${digit}0000s --root &
sleep 1
python plot_double_coincidences.py \
"${in_path}/EH3/hadded_ad1/out_ad1_${digit}????.root" \
${out_path}/EH3/double_coincidences_ad1/double_coincs_EH3_AD1_${digit}0000s --root &
sleep 1
python plot_double_coincidences.py \
"${in_path}/EH3/hadded_ad2/out_ad2_${digit}????.root" \
${out_path}/EH3/double_coincidences_ad2/double_coincs_EH3_AD2_${digit}0000s --root &
sleep 1
python plot_double_coincidences.py \
"${in_path}/EH3/hadded_ad3/out_ad3_${digit}????.root" \
${out_path}/EH3/double_coincidences_ad3/double_coincs_EH3_AD3_${digit}0000s --root &
sleep 1
python plot_double_coincidences.py \
"${in_path}/EH3/hadded_ad4/out_ad4_${digit}????.root" \
${out_path}/EH3/double_coincidences_ad4/double_coincs_EH3_AD4_${digit}0000s --root &
sleep 1
python plot_dr_vs_dt.py \
"${in_path}/EH1/hadded_ad1/out_ad1_${digit}????.root" \
${out_path}/EH1/dr_vs_dt_ad1/dr_vs_dt_EH1_AD1_${digit}0000s --root &
sleep 1
python plot_dr_vs_dt.py \
"${in_path}/EH1/hadded_ad2/out_ad2_${digit}????.root" \
${out_path}/EH1/dr_vs_dt_ad2/dr_vs_dt_EH1_AD2_${digit}0000s --root &
sleep 1
python plot_dr_vs_dt.py \
"${in_path}/EH2/hadded_ad1/out_ad1_${digit}????.root" \
${out_path}/EH2/dr_vs_dt_ad1/dr_vs_dt_EH2_AD1_${digit}0000s --root &
sleep 1
python plot_dr_vs_dt.py \
"${in_path}/EH2/hadded_ad2/out_ad2_${digit}????.root" \
${out_path}/EH2/dr_vs_dt_ad2/dr_vs_dt_EH2_AD2_${digit}0000s --root &
sleep 1
python plot_dr_vs_dt.py \
"${in_path}/EH3/hadded_ad1/out_ad1_${digit}????.root" \
${out_path}/EH3/dr_vs_dt_ad1/dr_vs_dt_EH3_AD1_${digit}0000s --root &
sleep 1
python plot_dr_vs_dt.py \
"${in_path}/EH3/hadded_ad2/out_ad2_${digit}????.root" \
${out_path}/EH3/dr_vs_dt_ad2/dr_vs_dt_EH3_AD2_${digit}0000s --root &
sleep 1
python plot_dr_vs_dt.py \
"${in_path}/EH3/hadded_ad3/out_ad3_${digit}????.root" \
${out_path}/EH3/dr_vs_dt_ad3/dr_vs_dt_EH3_AD3_${digit}0000s --root &
sleep 1
python plot_dr_vs_dt.py \
"${in_path}/EH3/hadded_ad4/out_ad4_${digit}????.root" \
${out_path}/EH3/dr_vs_dt_ad4/dr_vs_dt_EH3_AD4_${digit}0000s --root &
sleep 1
python plot_post_DT_cut.py \
"${in_path}/EH1/hadded_ad1/out_ad1_${digit}????.root" \
${out_path}/EH1/post_DT_cut_ad1/post_DT_cut_EH1_AD1_${digit}0000s --root &
sleep 1
python plot_post_DT_cut.py \
"${in_path}/EH1/hadded_ad2/out_ad2_${digit}????.root" \
${out_path}/EH1/post_DT_cut_ad2/post_DT_cut_EH1_AD2_${digit}0000s --root &
sleep 1
python plot_post_DT_cut.py \
"${in_path}/EH2/hadded_ad1/out_ad1_${digit}????.root" \
${out_path}/EH2/post_DT_cut_ad1/post_DT_cut_EH2_AD1_${digit}0000s --root &
sleep 1
python plot_post_DT_cut.py \
"${in_path}/EH2/hadded_ad2/out_ad2_${digit}????.root" \
${out_path}/EH2/post_DT_cut_ad2/post_DT_cut_EH2_AD2_${digit}0000s --root &
sleep 1
python plot_post_DT_cut.py \
"${in_path}/EH3/hadded_ad1/out_ad1_${digit}????.root" \
${out_path}/EH3/post_DT_cut_ad1/post_DT_cut_EH3_AD1_${digit}0000s --root &
sleep 1
python plot_post_DT_cut.py \
"${in_path}/EH3/hadded_ad2/out_ad2_${digit}????.root" \
${out_path}/EH3/post_DT_cut_ad2/post_DT_cut_EH3_AD2_${digit}0000s --root &
sleep 1
python plot_post_DT_cut.py \
"${in_path}/EH3/hadded_ad3/out_ad3_${digit}????.root" \
${out_path}/EH3/post_DT_cut_ad3/post_DT_cut_EH3_AD3_${digit}0000s --root &
sleep 1
python plot_post_DT_cut.py \
"${in_path}/EH3/hadded_ad4/out_ad4_${digit}????.root" \
${out_path}/EH3/post_DT_cut_ad4/post_DT_cut_EH3_AD4_${digit}0000s --root &
sleep 1
python plot_singles.py \
"${in_path}/EH1/singles_ad1/out_ad1_${digit}????.root" \
${out_path}/EH1/singles_ad1/singles_EH1_AD1_${digit}0000s --root &
sleep 1
python plot_singles.py \
"${in_path}/EH1/singles_ad2/out_ad2_${digit}????.root" \
${out_path}/EH1/singles_ad2/singles_EH1_AD2_${digit}0000s --root &
sleep 1
python plot_singles.py \
"${in_path}/EH2/singles_ad1/out_ad1_${digit}????.root" \
${out_path}/EH2/singles_ad1/singles_EH2_AD1_${digit}0000s --root &
sleep 1
python plot_singles.py \
"${in_path}/EH2/singles_ad2/out_ad2_${digit}????.root" \
${out_path}/EH2/singles_ad2/singles_EH2_AD2_${digit}0000s --root &
sleep 1
python plot_singles.py \
"${in_path}/EH3/singles_ad1/out_ad1_${digit}????.root" \
${out_path}/EH3/singles_ad1/singles_EH3_AD1_${digit}0000s --root &
sleep 1
python plot_singles.py \
"${in_path}/EH3/singles_ad2/out_ad2_${digit}????.root" \
${out_path}/EH3/singles_ad2/singles_EH3_AD2_${digit}0000s --root &
sleep 1
python plot_singles.py \
"${in_path}/EH3/singles_ad3/out_ad3_${digit}????.root" \
${out_path}/EH3/singles_ad3/singles_EH3_AD3_${digit}0000s --root &
sleep 1
python plot_singles.py \
"${in_path}/EH3/singles_ad4/out_ad4_${digit}????.root" \
${out_path}/EH3/singles_ad4/singles_EH3_AD4_${digit}0000s --root &
wait
