# Investigating correlations between PFS and mutational/demographical feautures based on Liu, et al. (Nature Communications, 2017)

**Author**: Leonidas Tsepenekas ([ltsepene@cs.umd.edu](mailto:ltsepene@cs.umd.edu))
**Date**: November 14, 2018

## Plan of investigation

Liu, et al. [1] performed discovery of mutational signatures analysis in chemotherapy resistant muscle-invasive bladder cancer for three cohorts. The first cohort involved pre-treatment mutations, the second post-treatment mutations and the last post-only mutations for the same set of patients. Their study also includes information about the progression-free survival (PFS) period for each patient. In this project, our goal is investigate interesting factors that may correlate with PFS. In the discussion that follows, we measure PFS as the days passed between cystectomy and the day the cancer reoccured. 

We performed the following set of experiments:
1. For each patient we computed the pre-treatment and post-only treatment signature profile, using the NMF approach of Liu et al., that we implemented for project 1. The goal of this experiment was to identify any relationship between pre and post-only signature profiles.
2. We tried to associate PFS with specific signatures in post-only treatment signature profiles. 
3. We tried to associate PFS with specific signatures in pre-treatment signature profiles. 
4. We tried to associate PFS with the age and sex of the patients in our cohort.
5. We tried to associate PFS with mutational load in all three cases considered. That is for, pre-treatment, post-treatment, and post-only treatment mutations.

## Data
Initially, the data we used included the mutational counts given to us for the first project. Those data include 30 samples for the pre-treatment and post-treatment case, but only 29 for the post-only case. We used those datasets to calculate signatures and exposures, as in project 1, but also mutational burdens. All the additional data required were available in Supplementary material of the Liu et al. submission website https://www.nature.com/articles/s41467-017-02320-7#Sec29 (Supplementary Data 1). To match demographic features and PFS information to each patient we used the patient ids provided in both the datasets of project 1 and the above supplementary table. Furthermore, we had to identify which patient was missing in the post-only case and make the necessary changes in our code in order to correctly handle that.

## Experiment 1 
For this experiment we tried to infer if there are any interesting correlations between pre-treatment and post-only treatment signature profiles. We chose to study the post-only scenario, since our NMF implementation from project 1 was able to fully reproduce the results for this case, compared to the post-treatment one. Our goal was to possibly identify if there exists a pre-treatment profile that strongly correlates with increased activity of the newly discovered cisplatin signature. 

To compute signatures and exposures we used our code from project 1, setting the effective number of signatures to the values used in Liu et al. The only difference is that this time we used more NMF runs, specifically 300 instead of 100, in order to get more accurate results. The following graph shows the exposures for both cases for all 29 patients of interest. Each bar corresponds to a patient (the ordering of the patients is the same in both graphs), and the bar plots show the activity of each signature.

![Experiment 1](exp1.PNG)

It is clear, even to the naked eye, that no correlation between pre and post-only exposures exists. Take for instance the pairs of patients (5,6) and (28,29) (the last two). Both patients in each pair have remarkably similart pre-treatment exposures, but quite different post-only exposures. Moreover, the activity of the cisplatin signatures doesn't seem to depend on the activity of either the APOBEC or the NER signature.

## Experiment 2
For this experiment we tried to infer correlations between post-only exposures and PFS. The exposures were calculated as mentioned earlier. We stratified patients based on their dominant signature, that is the signature with the highest exposure value among the three (APOBEC, AGING, Cisplatin). We also used cosine similarities in order to differentiate the discovered signatures. The following box plots show our results. The orange line indicates the median value and circles correspond to outliers. Also, the first box is for the APOBEC signature, the second for the AGING and the third for the cisplatin one.

![Experiment 2](exp2.PNG)


### Inferring the effective number of signatures

The figures below show the statistics we got. EV stands for expected variance, CC for cophenetic coefficients, RSS for residual sum of squares, RE for residual error and
ACS for average cosine similarity.

Pre-treatment:
![Pre-treatment](pre-k.jpg)

Post-treatment:
![Post-treatment](post-k.jpg)

Post-only-treatment:
![Post only](post-only-k.jpg)

For the pre-treatment case Liu et al. chose k=2. However, in our statistics only 2 out of the 4 metrics suggest that as a good choice. Namely, the cophenetic coefficients 
and the cosine similarity. Given that Liu et al. use the Brunet update method in their NMF, which relies strongly on cophenetic coefficients, k=2 may not be that bad after all.
For the post-treatment case the authors use k=4. We believe that our measurements agree with that. The average cosine similarity is maximized for that value and also all
other metrics take reasonable values there. Finally, in the post-only case, Liu et al. choose k=3. To our understanding, our results also agree with, namely for the above reasons.
Therefore, we can conclude that in some sense, that inferring the number of active signatures is reproducible.

### Signature Discovery

If we use the number of signatures used by the authors we get surprisingly similar results. Some minor differences are observed in the second cohort only. 

Liu's signatures:
![Liu's et al. results](Liu.jpg) 

Pre-treatment signatures we discovered:
![Pre-treatment](pre.jpg)

Post-treatment signatures we discovered:
![Post-treatment](post.jpg)

Post-only-treatment signatures we discovered:
![Post only](post-only.jpg)

The only differences are those observed in the second cohort. For that case, we are able to discover 2 out of the 4 signatures presented in the paper (the last two in
our figure). Our first signature looks like a combination of the of signatures 3 and 4 of the paper. As for our second signature, we can't see any similarity between that
and any of the signatures in the paper. That's mainly because of the many C > A mutations. For the other two cohorts, our results look suprisingly similar to that of the 
paper. The only difference lies in the cosine similarities with the cosmic collection. In some cases, we get the same Cosmic signature correspondings as the authors. But 
in other cases we get a combination of signatures instead of a single cosmic one. This may be the case because the authors may not have made explicit comparisons with combinations
of cosmic signatures for all of the discovered ones. Overall, reproducing the experiment given the same k as in the paper seems possible. Minor differences may result from
details in the use of NMF that the authors don't mention.

## References
1. Liu, et al. (2017) "Mutational patterns in chemotherapy resistant muscle-invasive bladder cancer." _Nature Communications_ **8**, Article number: 2193. [doi: 10.1038/s41467-017-02320-7](https://doi.org/10.1038/s41467-017-02320-7)
