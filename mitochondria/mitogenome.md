# Mitochondria analysis

```bash
muscle -in curryco_mt.fasta -out Pra_mytogenome_aln.fasta -maxiters 3

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
