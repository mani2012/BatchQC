---
title: "BatchQC Report"
date: "2017-11-27"
output: 
  html_vignette:
    toc: true
    toc_depth: 2
    template: batchQC.html
    self_contained: no
    lib_dir: libs
---


Summary
=======
## Confounding
### Number of samples in each Batch and Condition

------------------------------------------------------------------------
        &nbsp;          Batch 1   Batch 2   Batch 3   Batch 4   Batch 5 
---------------------- --------- --------- --------- --------- ---------
 **Condition Biopsy**      0         0         0         5         4    

 **Condition Cancer**     11        14         0         0        15    

 **Condition Normal**      0         4         4         0         0    
------------------------------------------------------------------------

### Measures of confounding between Batch and Condition

----------------------------------------------------------------------
            &nbsp;                Standardized Pearson     Cramer's V 
                                 Correlation Coefficient              
------------------------------- ------------------------- ------------
  **Confounding Coefficients             0.8871              0.7428   
 (0=no confounding, 1=complete                                        
        confounding)**                                                
----------------------------------------------------------------------

## Variation Analysis
### Variation explained by Batch and Condition
![](C:\Users\zhang\Dropbox\Work\EvanJohnson\refComBat_paper_092116\BatchQC\combat_batchqc_report_files/figure-html/unnamed-chunk-4-1.png)<!-- -->


----------------------------------------------------------
   &nbsp;      Full (Condition+Batch)   Condition   Batch 
------------- ------------------------ ----------- -------
  **Min.**             0.113              0.002     0.106 

 **1st Qu.**           21.72              11.55     11.72 

 **Median**            37.29              27.54     18.48 

  **Mean**             37.89              29.61     18.43 

 **3rd Qu.**           53.68              46.1      24.78 

  **Max.**             92.54              91.05     57.11 
----------------------------------------------------------

## P-value Analysis
### Distribution of Batch and Condition Effect p-values Across Genes

----------------------------------------------------------------------------------------------
         &nbsp;             Min.      1st Qu.     Median     Mean    3rd Qu.   Max.   Ps<0.05 
------------------------ ---------- ----------- ---------- -------- --------- ------ ---------
   **Batch P-values**     6.94e-08    0.05631     0.2664    0.3509    0.609     1     0.2374  

 **Condition P-values**      0       3.689e-06   0.001872   0.1227   0.1091     1      0.69   
----------------------------------------------------------------------------------------------

![](C:\Users\zhang\Dropbox\Work\EvanJohnson\refComBat_paper_092116\BatchQC\combat_batchqc_report_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

![](C:\Users\zhang\Dropbox\Work\EvanJohnson\refComBat_paper_092116\BatchQC\combat_batchqc_report_files/figure-html/unnamed-chunk-8-1.png)<!-- -->


Differential Expression
=======================
## Expression Plot
Boxplots for all values for each of the samples and are colored by batch membership.

![](C:\Users\zhang\Dropbox\Work\EvanJohnson\refComBat_paper_092116\BatchQC\combat_batchqc_report_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

## LIMMA

-------------------------------------------------------------------------------------------------------------------
     &nbsp;        Condition: Cancer (logFC)   Condition: Normal (logFC)   AveExpr     F      P.Value    adj.P.Val 
----------------- --------------------------- --------------------------- --------- ------- ----------- -----------
  **216005_at**             -3.338                      -2.759              4.289     173    5.03e-24    1.121e-19 

  **205207_at**             -4.886                      -4.524              5.344     144    3.258e-22   3.63e-18  

  **220287_at**              -2.8                       -1.773              6.081    121.7   1.36e-20    1.01e-16  

  **211585_at**             -2.121                      -0.7714             4.596    99.35   1.058e-18   5.893e-15 

  **220874_at**             -2.421                      -0.9097             5.369    94.09   3.293e-18   1.467e-14 

  **220244_at**             -2.547                      -1.097              6.467    91.46   5.93e-18    2.202e-14 

 **222329_x_at**            -3.514                      0.3567              7.12     88.35   1.205e-17   3.625e-14 

 **211143_x_at**             -2.02                      -1.315              6.411    88.02   1.302e-17   3.625e-14 

  **211565_at**             -3.882                      -1.114              5.865    87.37   1.514e-17   3.749e-14 

 **209057_x_at**            -3.438                      0.08605             8.162    86.09   2.048e-17   4.563e-14 
-------------------------------------------------------------------------------------------------------------------


Median Correlations
===================
This plot helps identify outlying samples.
![](C:\Users\zhang\Dropbox\Work\EvanJohnson\refComBat_paper_092116\BatchQC\combat_batchqc_report_files/figure-html/unnamed-chunk-13-1.png)<!-- -->


Heatmaps
========
## Heatmap
This is a heatmap of the given data matrix showing the batch effects and variations with different conditions.
![](C:\Users\zhang\Dropbox\Work\EvanJohnson\refComBat_paper_092116\BatchQC\combat_batchqc_report_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

## Sample Correlations
This is a heatmap of the correlation between samples.
![](C:\Users\zhang\Dropbox\Work\EvanJohnson\refComBat_paper_092116\BatchQC\combat_batchqc_report_files/figure-html/unnamed-chunk-16-1.png)<!-- -->


Circular Dendrogram
===================
This is a Circular Dendrogram of the given data matrix colored by batch to show the batch effects.
![](C:\Users\zhang\Dropbox\Work\EvanJohnson\refComBat_paper_092116\BatchQC\combat_batchqc_report_files/figure-html/unnamed-chunk-18-1.png)<!-- -->


PCA: Principal Component Analysis
=================================
## PCA
This is a plot of the top two principal components colored by batch to show the batch effects.
![](C:\Users\zhang\Dropbox\Work\EvanJohnson\refComBat_paper_092116\BatchQC\combat_batchqc_report_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

## Explained Variation

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  &nbsp;    Proportion of Variance (%)   Cumulative Proportion of   Percent Variation Explained by   Percent Variation Explained by   Condition Significance   Percent Variation Explained by   Batch Significance (p-value) 
                                               Variance (%)           Either Condition or Batch                Condition                    (p-value)                      Batch                                             
---------- ---------------------------- -------------------------- -------------------------------- -------------------------------- ------------------------ -------------------------------- ------------------------------
 **PC1**              38.03                       38.03                           73                              65.8                          0                           28.2                          0.01621            

 **PC2**              9.321                       47.35                          39.9                             20.6                       0.08866                        33.8                           0.0067            

 **PC3**              5.105                       52.45                           32                              5.3                         0.8358                        31.5                          0.00201            

 **PC4**              3.425                       55.88                          20.8                              17                         0.1887                        15.4                           0.6564            

 **PC5**              3.026                       58.91                          37.2                              27                        0.00047                        14.6                           0.1039            

 **PC6**              2.484                       61.39                          19.6                             16.7                       0.05088                        9.4                            0.7759            

 **PC7**              2.049                       63.44                          3.6                              2.1                         0.7844                        2.7                            0.9378            

 **PC8**              1.964                        65.4                           28                               6                         0.03409                        17.5                           0.0088            

 **PC9**              1.775                       67.18                          10.1                             8.3                         0.3409                        6.1                            0.9084            

 **PC10**             1.707                       68.88                           11                              4.5                         0.526                         8.6                            0.4662            

 **PC11**             1.536                       70.42                           7                               1.4                         0.6694                        5.4                            0.5645            

 **PC12**             1.419                       71.84                          2.3                              0.3                         0.8115                        1.4                            0.9071            

 **PC13**             1.359                        73.2                          6.7                              0.8                         0.6956                        5.3                            0.5378            

 **PC14**             1.208                       74.41                          4.8                              3.5                         0.452                         1.7                            0.955             

 **PC15**              1.17                       75.58                          8.2                              0.7                         0.9527                        8.1                            0.4058            

 **PC16**             1.083                       76.66                          2.7                              0.1                         0.9623                        2.5                            0.8563            

 **PC17**             1.051                       77.71                          4.7                              0.7                         0.7809                        3.7                            0.7228            

 **PC18**             1.012                       78.72                          2.4                              0.6                         0.7658                        1.4                            0.9172            

 **PC19**             0.9901                      79.71                          5.3                              0.7                         0.868                         4.7                            0.6673            

 **PC20**             0.9734                      80.69                          5.9                              3.3                         0.4364                        2.7                            0.8505            

 **PC21**             0.9076                      81.59                          1.5                              0.3                         0.8675                         1                             0.9607            

 **PC22**             0.8812                      82.47                          1.9                              0.5                         0.8619                        1.3                            0.9499            

 **PC23**             0.849                       83.32                          1.4                              0.6                         0.8273                        0.6                            0.983             

 **PC24**             0.8135                      84.14                          8.3                              0.5                         0.9526                        8.1                            0.385             

 **PC25**             0.7909                      84.93                           2                               0.1                         0.8976                        1.6                            0.9111            

 **PC26**             0.7806                      85.71                          4.2                              0.6                         0.7536                        3.1                            0.759             

 **PC27**             0.7475                      86.46                          1.4                              0.3                         0.8143                        0.6                            0.9628            

 **PC28**             0.6878                      87.14                           3                               0.6                         0.9934                        2.9                            0.8701            

 **PC29**             0.6675                      87.81                          1.9                              0.5                         0.7774                        0.9                            0.9471            

 **PC30**             0.6507                      88.46                          4.5                              0.1                         0.9128                        4.1                            0.6864            

 **PC31**             0.643                       89.11                          1.6                              0.1                         0.9022                        1.2                            0.9425            

 **PC32**             0.6309                      89.74                          5.7                              0.4                         0.806                         4.9                            0.5892            

 **PC33**             0.6068                      90.34                          2.6                              1.4                         0.7175                        1.2                            0.9619            

 **PC34**             0.5764                      90.92                           4                               0.5                         0.6021                        2.1                            0.7637            

 **PC35**             0.5689                      91.49                          1.1                              0.1                         0.9781                         1                             0.9713            

 **PC36**             0.5522                      92.04                          2.9                              0.3                         0.7218                        1.6                            0.8568            

 **PC37**             0.5342                      92.57                          2.3                               0                          0.9962                        2.3                            0.8768            

 **PC38**             0.5238                       93.1                          4.9                              0.4                         0.8891                        4.5                            0.6643            

 **PC39**             0.5102                      93.61                          4.4                              0.5                         0.9883                        4.4                            0.7284            

 **PC40**             0.4961                       94.1                          1.7                               0                          0.9653                        1.5                            0.9329            

 **PC41**             0.4955                       94.6                          3.2                              0.1                         0.8694                        2.7                            0.8005            

 **PC42**             0.4734                      95.07                          0.4                              0.3                         0.9534                        0.2                            0.9996            

 **PC43**             0.4684                      95.54                          0.8                              0.4                         0.9932                        0.8                            0.9941            

 **PC44**             0.4527                      95.99                          1.7                              0.9                         0.8072                        0.9                            0.9805            

 **PC45**             0.4501                      96.44                          10.1                             0.4                         0.9276                        9.8                            0.2679            

 **PC46**             0.4219                      96.87                          0.8                               0                          0.9072                        0.4                            0.9828            

 **PC47**             0.4128                      97.28                           3                                0                          0.8847                        2.5                            0.8238            

 **PC48**             0.3939                      97.67                          0.8                              0.1                         0.9478                        0.5                            0.9873            

 **PC49**             0.3762                      98.05                          1.6                              0.1                         0.9142                        1.2                            0.9412            

 **PC50**             0.3544                       98.4                          1.6                              0.1                         0.9627                        1.4                            0.9419            

 **PC51**             0.3449                      98.75                          6.6                              0.2                         0.9572                        6.4                            0.5001            

 **PC52**             0.3319                      99.08                           7                               0.1                         0.894                         6.6                            0.4503            

 **PC53**             0.3011                      99.38                           2                                0                          0.9606                        1.8                            0.9066            

 **PC54**             0.2663                      99.65                          0.9                              0.1                         0.8956                        0.5                            0.9818            

 **PC55**             0.2365                      99.88                          77.1                             0.9                           0                           49.3                             0               

 **PC56**             0.1152                       100                           96.2                             3.3                           0                            61                              0               

 **PC57**            7.38e-29                      100                           45.9                             39.7                        7e-05                         20.7                           0.2347            
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


Shape
=====
This is a heatmap plot showing the variation of gene expression mean, variance, skewness and kurtosis between samples grouped by batch to see the batch effects variation
![](C:\Users\zhang\Dropbox\Work\EvanJohnson\refComBat_paper_092116\BatchQC\combat_batchqc_report_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

```
## Note: Sample-wise p-value is calculated for the variation across samples on the measure across genes. Gene-wise p-value is calculated for the variation of each gene between batches on the measure across each batch. If the data is quantum normalized, then the Sample-wise measure across genes is same for all samples and Gene-wise p-value is a good measure.
```


Combat Plots
============
This is a plot showing whether parametric or non-parameteric prior is appropriate for this data. It also shows the Kolmogorov-Smirnov test comparing the parametric and non-parameteric prior distribution.

```
## Found 5 batches
## Adjusting for 2 covariate(s) or covariate level(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
```

![](C:\Users\zhang\Dropbox\Work\EvanJohnson\refComBat_paper_092116\BatchQC\combat_batchqc_report_files/figure-html/unnamed-chunk-25-1.png)<!-- -->![](C:\Users\zhang\Dropbox\Work\EvanJohnson\refComBat_paper_092116\BatchQC\combat_batchqc_report_files/figure-html/unnamed-chunk-25-2.png)<!-- -->![](C:\Users\zhang\Dropbox\Work\EvanJohnson\refComBat_paper_092116\BatchQC\combat_batchqc_report_files/figure-html/unnamed-chunk-25-3.png)<!-- -->![](C:\Users\zhang\Dropbox\Work\EvanJohnson\refComBat_paper_092116\BatchQC\combat_batchqc_report_files/figure-html/unnamed-chunk-25-4.png)<!-- -->

```
## Batch mean distribution across genes: Normal vs Empirical distribution
## Two-sided Kolmogorov-Smirnov test
## Selected Batch: 1
## Statistic D = 0.01064
## p-value = 0.01286
## 
## 
## Batch Variance distribution across genes: Inverse Gamma vs Empirical distribution
## Two-sided Kolmogorov-Smirnov test
## Selected Batch: 1
## Statistic D = 0.08217
## p-value = 0Note: The non-parametric version of ComBat takes much longer time to run and we recommend it only when the shape of the non-parametric curve widely differs such as a bimodal or highly skewed distribution. Otherwise, the difference in batch adjustment is very negligible and parametric version is recommended even if p-value of KS test above is significant.
```


SVA
===
## Summary

```
## Number of Surrogate Variables found in the given data: 10
```
