mkdir intensities
for j in $(seq 1 4); do mv B004_CL_segmentation/quantifications/reg00"$j"_statistics_growth_3_comp.csv intensities/B004_CL_reg_"$j".csv ; done
for j in $(seq 1 4); do mv B004_SB_segmentation/quantifications/reg00"$j"_statistics_growth_3_comp.csv intensities/B004_SB_reg_"$j".csv ; done
for i in 10 11; do for j in $(seq 1 4); do mv B0"$i"A_segmentation/quantifications/reg00"$j"_statistics_growth_3_comp.csv intensities/B0"$i"A_reg_"$j".csv ; done; done
for i in 10 11; do for j in $(seq 1 4); do mv B0"$i"B_segmentation/quantifications/reg00"$j"_statistics_growth_3_comp.csv intensities/B0"$i"B_reg_"$j".csv ; done; done

#sed -i 's/PDPN/Podoplanin/' intensities/B005*
#sed -i 's/Synapt/Synapto/' intensities/B005*
sed -i 's/aDefensin5/aDef5/' intensities/B01*
sed -i 's/Synaptophysin/Synapto/' intensities/B01*