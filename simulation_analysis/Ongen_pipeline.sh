#ONGEN PIPELINE
sbatch -c 1 -N 1 -t 0-01:00 -p short --mem=2G -o out/run_Ongen_eqtl.out --wrap="Rscript Nov_Ongen_1_runGWAS.R"

#for every tissue, do genome-wide eqtl analysis. Get list of eqtls that are significant and set of null variants too. 
for sim in {1..1000}
do
sbatch -c 1 -N 1 -t 0-00:15 -p short --mem=2G -o out/run_Ongen_eqtl_${sim}.out --wrap="Rscript Nov_Ongen_1_runeQTL.R ${sim}"
done

#GWAS-eqtl colocalization
samp=0 #step 1 (b 3-12 contigent on b 1-2)
ve=0.1
for sim in {1..1000}
do
for b in {1..2}
do
sbatch -c 1 -N 1 -t 0-10:00 -p short --mem=2G -o out/run_OngenR_${sim}_${b}.out --wrap="Rscript Nov_Ongen_2_GWAS_eQTL_sharing_prob.R ${sim}_${b}_${samp}_${ve}" #coloc of GWAS and eqtl
done
done

#Tissue-tissue colocalization
ve=0.1
for sim in {1..1000}
do
for b in {3..12} #tissues 3-12 (0-9); Rscript iterates through all tissues
do
for samp in {1..6}
do
sbatch -c 1 -N 1 -t 0-0${samp}:00 -p short --mem=2G -o out/run_OngenR_${sim}_${b}_${samp}_real.out --wrap="Rscript Nov_Ongen_2_GWAS_eQTL_sharing_prob_v2_rerun.R ${sim}_${b}_${samp}_${ve}" #coloc of GWAS and eqtl
sbatch -c 1 -N 1 -t 0-0${samp}:00 -p short --mem=2G -o out/run_OngenR_${sim}_${b}_${samp}_null.out --wrap="Rscript Nov_Ongen_2_GWAS_eQTL_sharing_prob_v2_rerunnull.R ${sim}_${b}_${samp}_${ve}"
#need separate output files
done
done
done

#Compute sharing probability of tissue eqtls
for sim in {1..100}
do
sbatch -c 1 -N 1 -t 0-00:05 -p short --mem=2G -o out/run_OngenR_last_${sim}.out --wrap="Rscript Nov_Ongen_4_eQTL_eQTL_sharing_prob.R $sim"
done

#Assess power and Type 1 error
Rscript Nov_Ongen_analyze.R

