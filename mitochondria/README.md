# Mitochondria analysis

The _P. ramorum_ sequences were used to retrieved 295 mitochondrial genomes using a guided assembly, 160 for EU1 and 135 for NA1. The mitochondrial genomes were assembled and aligned with [clustal omega](http://www.clustal.org/omega/). The alignments were used to construct two phylogenetic trees, one using the Neighbor-joining method, and the second to reconstruct a maximum likelihood tree. 

##### Command line for using clustalo at the CQLS
```bash
$ clustalo -i Pram_mitogenome.fasta -o Pram_mitogenome_aln.fas --iter=3 --threads=20
# SGE submission EU1
SGE_Batch -c 'clustalo -i Pram_mitogenome_EU1.fasta -o Pram_mitogenome_EU1-aln.fas --iter=10 --threads=12 -v' -P 12 -q bpp -r clustal-EU1mg
# SGE submission NA1
SGE_Batch -c 'clustalo -i Pram_mitogenome_NA1.fasta -o Pram_mitogenome_NA1-aln.fas --iter=10 --threads=12 -v' -P 12 -q bpp -r clustal-NA1mg

```
> The genome lenght of the mitochondrial genome was 39,739 bp with a mean lenght of 39,317, 96% of identical sites and a pairwise identity of 99.%. The GC content rounds about 22%.
> Results: EU1 failed, running again on 30th. NA1 alignment done.

## Temporal signal test:
The Neighbor-joining tree was used to calculate the temporal signal using the program [TemPest](https://beast.community/tempest). 
The dated tips showed some temporal signal using the correlation function **Correlation Coefficient (CC):** 0.582 fitting the best root, and 0.527 without it. The CC for residual mean squared was 0.527, and was reduced to 0.309 while fitting the best root.\
**Correlation Function results:**
| Dated tips | no-root | best-fitting root |
|--|--|--|
| date range | 18     | 18 |
| slope(rate)| 2.004e-5| 2.420e-5 |
|X-interception | 1988 | 1965 |
|Correlation Coefficient | 0.5273| 0.5828|
|R Squared  | 0.2781 | 0.3397 |
|Residual Mean Square| 5.752e-8| 6.2831 |

![tempest_results](https://github.com/ricardoi/genedis_networks/blob/main/mitochondria/figures/Pram_295_tempest.png)

### Results for Lineages NA1 & EU1
I used a maximum likelihood tree to calculate the temporal signal using the program [TemPest](https://beast.community/tempest) for both lineages, NA1 and EU1. The results indicated that there is no temporal signal to calculate a molecular clock. 

**Correlation Function results:**

| Dated tips | no-root | best-fitting root |
|--|--|--|
|            | NA1     |                   |
| date range | 4     | 4 |
| slope(rate)| 1.51e-5| 9.846e-5 |
|X-interception | 1988 | 1906 |
|Correlation Coefficient | 7.806e-2 | 5.867e-2|
|R Squared  | 6.094e-3 | 3.443e-3 |
|Residual Mean Square| 9.082e-8| 6.772e-8 |
|            | EU1     |                   |
| date range | 3     | 3 |
| slope(rate)| 2.95e-5| 7.132e-5 |
|X-interception | 2005 | 1972 |
|Correlation Coefficient | 0.111| 3.888e-2|
|R Squared  | 1.239e-2 | 1.511e-3 |
|Residual Mean Square| 5.260e-8| 2.549e-8 |

![tempest_results2](https://github.com/ricardoi/genedis_networks/blob/main/mitochondria/figures/NA1-EU1_tempest_results.png)

##### Prior selection 
Important reading to set [prior selection](https://taming-the-beast.org/tutorials/Prior-selection/) in beast.


# Maximum Likelihood
IQTREE maximum likelihood, standard command `iqtree2 -s sequence_file.fas -m TEST -B 1000` testing for evolutionary model, with 1000 bootstraps.
```bash 
SGE_Batch -c 'iqtree2 -s Pram_mitogenome_NA1-aln.fas -m TEST -B 1000 -T 8 -v' -P 8 -q bpp -r IQTR-NA1
```

## Programs Troubleshooting and Installations
```bash
$ muscle -in curryco_mt.fasta -out Pra_mytogenome_aln.fasta -maxiters 3

MUSCLE v3.8.31 by Robert C. Edgar

http://www.drive5.com/muscle
This software is donated to the public domain.
Please cite: Edgar, R.C. Nucleic Acids Res 32(5), 1792-97.

curryco_mt 283 seqs, max length 39394, avg  length 39317
00:00:15  70 MB(-887%)  Iter   1  100.00%  K-mer dist pass 1
00:00:15  70 MB(-887%)  Iter   1  100.00%  K-mer dist pass 2
01:44:48  8206 MB(-103317%)  Iter   1  100.00%  Align node      
01:44:48  8212 MB(-103383%)  Iter   1  100.00%  Root alignment
03:30:01  8212 MB(-103386%)  Iter   2  100.00%  Refine tree   
03:30:01  8212 MB(-103386%)  Iter   2  100.00%  Root alignment
03:30:01  8212 MB(-103386%)  Iter   2  100.00%  Root alignment

*** ERROR ***  MSA::GetLetter(0/1, 0/904)='???'/4294967295
12595.935u 7.760s 3:30:08.71 99.9%	0+0k 23816+0io 1pf+0w
```

### Installing BEAGLE and BEAST

```bash 
cmake -DCMAKE_INSTALL_PREFIX:PATH=$HOME -DBUILD_OPENCL=OFF -DCMAKE_CC_COMPILER=/usr/bin/cuda-gcc -DCMAKE_CXX_COMPILER=/usr/bin/cuda-g++ ..
```
Installed in `$HOME`, and running from `/bin`, GPU capabilities ON.

