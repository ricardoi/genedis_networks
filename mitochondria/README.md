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

The Neighbor-joining tree was used to calculate the temporal signal using the program TemPest.




## Troubleshooting programs and installing them
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

*** ERROR ***  MSA::GetLetter(0/1, 0/904)='�'/4294967295
12595.935u 7.760s 3:30:08.71 99.9%	0+0k 23816+0io 1pf+0w
```

### Installing BEAGLE and BEAST

```bash 
cmake -DCMAKE_INSTALL_PREFIX:PATH=$HOME -DBUILD_OPENCL=OFF -DCMAKE_CC_COMPILER=/usr/bin/cuda-gcc -DCMAKE_CXX_COMPILER=/usr/bin/cuda-g++ ..
```
Installed in `$HOME`, and running from `/bin`, GPU capabilities ON.

