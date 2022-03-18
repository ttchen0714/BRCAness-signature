# BRCAness signature Manual
## 1. Within-sample relative expression orderings (REO)
Within-sample REO was defined as a binary relation between gene i and gene j, where Gi > Gj or Gi < Gj. Gi and Gj denote the expression values of gene i and gene j, respectively. A set of REOs is highly stable in OvCa with BRCA1/2 mutation and reversed in patients with wild type BRCA1/2, which reflects the feature of BRCAness.<br> 
![](https://github.com/ttchen0714/BRCAness-signature/Image/f1.tif)<br>
## 2. BRCAness signature (2-GPS)
The workflow of developing the REO-based BRCAness signature: <br>
(1) we constructed a discrete profile based on the REOs of gene pairs;<br>
(2) we defined the reversed gene pair with a stable REO pattern (Gi > Gj) in more than 99% of mutation samples, but with a significantly reversed REO pattern (Gi < Gj) in wild-type samples tested by Fisher’s test with p < 0.05; <br>
![](https://github.com/ttchen0714/BRCAness-signature/Image/f2.tif)<br>
(3) a greedy algorithm was used to generate a set of gene pairs with high coverage in wild-type samples with iteration conditions. First, the gene pair with the highest coverage was chosen as the seed. Second, a candidate gene pair was
taken into account that increases the most coverage and more than 1%.<br>
Finally, using the rank-based algorithm, we discovered a BRCAness signature consisting of 2 gene pairs, including CXCL1 > SV2A and LY9 > CHRNB3, using gene expression profiles of ovarian cancer from the Cancer Genome Atlas.<br>
![](https://github.com/ttchen0714/BRCAness-signature/Image/f3.tif)<br>
## 3. Application
```
RDtij=Rti-Rtj;
```
Where RDtij reprents the rank difference between gene i and gene j in sample t, Rti reprents the rank of gene i in sample t and Rtj reprents the rank of gene j in sample t.<br>
```
BRCAness score(t)=count(RDtij>0)/N.  
```
For the OvCa sample t, the **BRCAness score** was calculated as the proportion of gene pairs (Gi > Gj) 
with values of rank differences more than 0 among the total gene pairs of the BRCAness signature (N=2). 
According to a determination rule, samples were divided into **"BRCAness group"** (BRCAness score = 1) or **"non-BRCAness group"** (BRCAness score < 1). <br>

