References

- the orginal DESeq publication by Anders and Huber

  [Differential expression analysis for sequence count data](http://genomebiology.com/2010/11/10/R106)

- DESeq documentation on Biocoductor website

  [DESeq](http://bioconductor.org/packages/release/bioc/html/DESeq.html)

Note that in the future we are likely to transition to
[DESeq2](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html)


### Example usage:

```
R --vanilla --args --table=example.tab  --header="YES" --rownamesCol=1 --treatCountCols=2,3 --contrCountCols=4,5 --outfiBaseTab=example --outfiBasePNG=example --filterNoPolya "YES" < DESEQ/deseq.r 
```
The input table
[example.tab](https://github.com/getopt/DESEQ/blob/master/Data/example.tab) and
output files are ind `Data/` directory of this repository. The head of input
[example.tab](https://github.com/getopt/DESEQ/blob/master/Data/example.tab)
file: 

```
$ head example.tab

Gene          cont_rep1   cont_rep2    treat_rep1  treat_rep2 
128up           1742       1675           1603        312
14-3-3epsilon   5164       5478           5885        940
14-3-3zeta      250        575            481         129
140up           223        282            219         54
18w             27         25             65          13
26-29-p         9634       9026           10666       1750
2mit            0          0              4           0
312             123        264            230         66
4EHP            52         50             70          10
```

In the above example, two text ouput files are generated:

1. [example.deseq.Results_all.tab.SKIPPED_IDS](https://github.com/getopt/DESEQ/blob/master/Data/example.deseq.Results_all.tab.SKIPPED_IDS) 
    file with gene names that were sckipped from the analysis

```
$ head example.deseq.Results_all.tab.SKIPPED_IDS 

mtRNApol
5.8SrRNA
LSU-rRNA_Dme|Protostomia
LSU-rRNA_Hsa|Metazoa
SSU-rRNA_Dme|Protostomia
SSU-rRNA_Hsa|Metazoa
l(2)not
l(3)neo18
l(3)neo38
l(3)neo43
```

2. [example.deseq.Results_all.tab](https://github.com/getopt/DESEQ/blob/master/Data/example.deseq.Results_all.tab)
   this is a standard DESeq ouput table. Note how rows with `Inf` and `NA`
   values are skipped using `cat` and `grep` for viewing the table

```
$ cat example.deseq.Results_all.tab | grep -v -P '\tInf\t' | grep -v -P '\tNA\t' | head

id      baseMean            baseMeanA           baseMeanB           foldChange          log2FoldChange      pval                padj
CG12655 5.52461827927325    0.310844928378475   10.738391630168     34.5458157743957    5.11043907533905    0.209850294873236   1
MtnC    16.4025209776703    0.932534785135426   31.8725071702052    34.1783573956188    5.09501115795936    0.475634174481844   1
CG4630  5.0114623846732     0.310844928378475   9.71207984096792    31.2441315727139    4.96551333599042    0.246396570329004   1
Amy-d   21.3660182286171    1.75309866697042    40.9789377902638    23.3751462837404    4.54690348842952    0.321955305030729   1
Sdic4   3.638867490771      0.310844928378475   6.96689005316353    22.41275123743      4.48624785007082    0.529797611756606   1
CG17192 3.57033394671122    0.310844928378475   6.82982296504396    21.9718011828978    4.45758123734407    0.584648893407203   1
CG5778  3.57033394671122    0.310844928378475   6.82982296504396    21.9718011828978    4.45758123734407    0.584648893407203   1
CG4363  18.7099759644333    1.75309866697042    35.6668532618962    20.3450347284407    4.34660483794369    0.333426563456943   1
CG30108 6.50521297073474    0.621689856756951   12.3887360847125    19.9275184403658    4.31669015846702    0.237003124886577   1
```

3. [outfi.deseq.dispersion.png](https://github.com/getopt/DESEQ/blob/master/Data/example.deseq.dispersion.png)
    plot of disperition estimates 

4. [outfi.deseq.diffexpres.png](https://github.com/getopt/DESEQ/blob/master/Data/example.deseq.dispersion.png)
   plot of differential expression estimates



### Summary usage:
```
R --vanilla --args --table=<file.tab>
                   --header=<YES|NO>
                   --rownamesCol=<N>
                   --fitType=<paramteric|local>
                   --method=<pooled|blind>
                   --sharingMode=<maximum|fit-only>
                   --treatCountCols=<1-based column numbers, csv>
                   --contrCountCols=<1-based column numbers, csv>
                   --outfiBaseTAB=<path/nameBase>
                   --outfiBasePDF=<path/nameBase> 
                   --filterNoPolya
```

### Description of arguments:
```
   --table            :   file with first column as row names and two columns
                      :   providing counts (values will be rounded to integers)

   --header           :   "YES" or "NO" (default "NO"; header isn't important
                          as such, only important for correct reading of the table)
   --rownamesCol      :   column number for the column that provides rownames


  --treatCountCols    :   <1-based column numbers, csv>
  --contrCountCols    :   <1-based column numbers, csv>


   --method           :   "pooled" or "blind" method of fitting dispersion model. 
                      :    Use "pooled" if you have replicates, "blind" if you 
                      :    do not have replicates. Default "pooled"

   --sharingMode      :   "maximum" or "fit-only" mode of fitting dispersion 
                      :    model. Use "maximum" if you have replicates, 
                      :    "fit-only" if you do not have replicates. Default: "maximum"

   --fitType          :   "parametric" or "local" mode of fitting dispersion model.
                      :   "parametric" (default) is recommended. But it doesn't 
                      :    work sometimes for unknown reasons, suggested alternative 
                      :   "local". Default: "parametric" 


   --outfiBaseTAB     :   path/nameBase of for the output files: 
                      :                                ${base}.deseq.Results_all.tab
   --outfiBasePNG     :   path/nameBase of for output PNGs. Two will be written: 
                                                       ${base}.deseq.dispersion.png 
                                                       ${base}.deseq.diffexpres.png 


   --filterNoPolya    :  set to "YES" if to remove certain gene names from the
                         table (those that are matched by "toFilter") by reg exp 
                         matching to gene names: 
                                                    "tRNA" "rRNA" "\\)n"
                                                    "4.5S" "5S" "7S" "_rich",
                                                    "\\|U\\d+" "snoRNA" "mir-"
```
*NOTE: more filtering is currently commented out inside of the script: rows
that have 0-value in all columns are not included into the analysis. Use of these
filtering steps will  results with `Inf` and `NA` values*

### Acknowledgments 
Thanks to **Junho Hur** for providing example.tab data file
