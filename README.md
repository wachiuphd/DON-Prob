# DON-Prob
Probabilistic analysis of DON.
HDMI is derived from 4 uncertainty distributions.
1. Bayesian benchmark dose (BBMD) distribution
2. Interspecies toxicokinetic (TK) distribution based on animal TK studies and human TK simulation
3. Default interspecies toxicodynamic (TD) distribution
4. Combination of intraspecies TK and TD distributions

BEMI in urine is derived from 4 uncertainty distributions.
1. BBMD distribution
2. Interspecies TK distribution based on animal TK studies and human TK simulation
3. Default interspecies TD distribution
4. Combination of intraspecies TK and TD distributions

BEMI in blood is derived from 3 uncertainty distributions.
1. BBMD distribution adjusted to average concentration in blood
2. Default interspecies TD distribution
3. Intraspecies TD distribution

## BBMD results

input files:
- iverson-f-et-al-1995-female-bmds.csv
- iverson-f-et-al-1995-male-bmds.csv

original files in repository:
BBMD plot

## Interspecies animal TK extrapolaion

script: 
  interspecies TK plot

output: 
- Supplementary Figure S1

## Human TK modeling results

input files from InputData folder:
- AUC_dose_samples.csv
- AUC_Utot_samples.csv

original files in repository:
DON-PK-model

## Human TD variability using LCL results

input files from InputData folder:
-DON-EC10_pop_samples.csv

original files in repository:
DON-LCL

## Combining datasets probabilistically

script:
- HDMI calculations.R
- BE calculation.R

output:
- HDMI.samples.csv
- HDMI.quantiles.csv
- Urine BE.samples.
- BE.urine.quantiles.csv
- Blood BE.samples.csv
- BE.blood.quantiles.csv
  
## Risk characterization 
Incidence of effects of the median population and fraction of change of decreased body weight based on different sensitive individuals are derived based on HDMI and BEMI functions. The incidence of effects of the median population prediction are compared with human biomonitoring data from Martins et al. (2019) and Wang et al. (2019).

input files:
- iverson-f-et-al-1995-female-bmds.csv - original files in repository: BBMD plot
- iverson-f-et-al-1995-male-bmds.csv - original files in repository: BBMD plot
- HDMI.samples.csv
- Urine BE.samples.csv
- Blood BE.samples.csv
- don exposure.csv - based on Martins et al. (2019) human biomonitoring exposure data, repository in "InputData" folder

script:
  Incidence of effect.Rmd
  
output:
- HDMI fraction of response.pdf - Incidence of effect and fraction of change of decreased body weight on different in sensitive populations based on HDMI function
- Urinary BE fraction of response.pdf - Incidence of effect and fraction of change of decreased body weight on different in sensitive populations based on urinary BEMI function
- Blood BE fraction of response.pdf - Incidence of effect and fraction of change of decreased body weight on different in sensitive populations based on blood BEMI function
- Log scale Population Incidence of 5 percent Decrease in Body Weight.pdf (Manuscript Figure 6)
- Percent change in BW for fractions of population.pdf (Supplementary Figure S8)
  
## Risk characterization
Individual Margin of Exposure is derived for sensitive population and random individuals for two study populations. Fraction of change of decreased body weight based on BEMI functions is calculated for different sensitive populations. Population exceeding TDI is calculated using EFSA (2017) TDI, WHO/IPCS (2018), the derived probabilitic TDI in this project.

input files:
- iverson-f-et-al-1995-female-bmds.csv - original files in repository: BBMD plot
- iverson-f-et-al-1995-male-bmds.csv - original files in repository: BBMD plot
- Urine BE.samples.csv
- Blood BE.samples.csv
- don exposure.csv - based on Martins et al. (2019) human biomonitoring exposure data, repository in "InputData" folder

script:
  Biomonitoring exposure.Rmd

output: 
- Martins et al (2019) population MOE quantile.csv
- Wang et al (2019) population MOE quantile.csv
- Percent Martins et al (2019) population UE fraction of response.csv
- Percent Wang et al (2019) population UE fraction of response.csv
- Martins et al (2019) Density plot_Fraction of response.pdf
- Quantiles of Martins et al (2019) population UE fraction of response.csv
- Quantiles of Wang et al (2019) population UE fraction of response.csv
- Boxplot Martins et al (2019) population UE Fraction of response.pdf
- Boxplot Wang et al (2019) population UE Fraction of response.pdf
- percent of population over TDI.csv