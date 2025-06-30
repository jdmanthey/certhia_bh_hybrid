# interactive job

# in windows directory

# combine the output for different analyses into a single file each
# first add a header for each file 
grep 'pop1' Ca_0001_Tg_2__100000001__100100000__stats.txt > ../window_heterozygosity.txt
grep 'pop1' Ca_0001_Tg_2__100000001__100100000__stats.txt > ../window_fst.txt
grep 'pop1' Ca_0001_Tg_2__100000001__100100000__stats.txt > ../window_dxy.txt
grep 'pop1' Ca_0001_Tg_2__100000001__100100000__stats.txt > ../window_titv.txt



# add the relevant stats to each file 
for i in $( ls *stats.txt ); do grep 'heterozygosity' $i >> ../window_heterozygosity.txt; done
for i in $( ls *stats.txt ); do grep 'Fst' $i >> ../window_fst.txt; done
for i in $( ls *stats.txt ); do grep 'Dxy' $i >> ../window_dxy.txt; done
for i in $( ls *stats.txt ); do grep 'titv' $i >> ../window_titv.txt; done




