# Reproducing Liu, et al. (_Nature Communications_, 2017)

**Author**: Leonidas Tsepenekas ([mdml@cs.umd.edu](mailto:mdml@cs.umd.edu))
**Date**: November 3, 2018

## Plan of investigation

Liu, et al. [1] performed discovery of mutational signatures analysis in in chemotherapy resistant muscle-invasive bladder cancer in three cohorts. 
The first cohort involved pre-treatment mutations, the second post-treatment mutations and the last post-only mutations. 

We reproduced this analysis, following the details described in the paper, with two goals in mind:
1. Infer the number of effective signatures for each cohort using the approach and metrics described in the paper. Vital details were missing, such as how they combine the 
results of their studied metrics when some of them disagree, and so we had to interpret the measurements on our own.
2. Given the above, implement the procedure described in the paper in order to actually determine the mutational signatures in each case.

### Infer number of effective signatures
The optimal rank (number of mutational signatures) was inferred after manually examining cophenetic coefficients, residuals, and residual sum of
squares for 50 NMF runs at ranks 2–8, as well as comparing discovered signatures to previously discovered signatures using a cosine similarity measure. High cophenetic coefficients, low residuals, low residual sum of squares, and high cosine
similarity to previous signatures were preferred. In addition to the above metrics used by Liu, et al., we also consider the expected variance. Moreover, since the paper 
frequently uses the average of previously discovered signatures in order to compute cosine similarity, we extended the cosmic signatures collection to a 30 * 30 collection of
all possible combinations of signatures, taking the average of the two each time. This new collection obviously contains the original cosmic signatures. We performed all
our comparisons based on that. Though the obtained results weren't much different from simply comparing with the original collection, we observed that this approach gave 
slightly tighter results for figuring out the effective number of signatures. Finally, statistics for residual sum of squares, cophenetic coefficients and expected variance 
were collected using the build in function estimate_rank() of nimfa. For the rest we manually collected residual errors and cosine similarities. For the residual error,
we computed the average error for each possible k and for cosine similarity we found for each discovered signature the best fitting choice for it (with respect to the known
combinations of cosmic signature) and then too the average of the highest cosine similarity of all discovered signatures.

#### Signature Discovery
As described in the paper, we performed 200 independent NMF runs for the inferred rank and
chose the resulting mutational signatures from the NMF run with the minimum residual error. We used the nimfa package for python in order to execute NMF. For the discovered
signatures we also computed their cosine similarity with the extended 30*30 cosmic collection. 

## Results, conclusions, and caveats

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