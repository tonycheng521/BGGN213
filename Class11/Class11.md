Protein Structure: bio3d package
================

``` r
summary = read.csv("Data Export Summary.csv")
summary
```

    ##   Experimental.Method Proteins Nucleic.Acids Protein.NA.Complex Other  Total
    ## 1               X-Ray   131463          2060               6768     8 140299
    ## 2                 NMR    11241          1304                262     8  12815
    ## 3 Electron Microscopy     2925            32               1004     0   3961
    ## 4               Other      280             4                  6    13    303
    ## 5        Multi Method      144             5                  2     1    152

``` r
total_protein = sum(summary["Proteins"])
x_ray_percent = summary[1,2] / total_protein *100
electron_micro_percent = summary[3,2] / total_protein *100
```

x-ray percentage is 90.0104756 %

``` r
sum(summary$Total)
```

    ## [1] 157530

``` r
ans = summary$Total/sum(summary$Total) * 100
names(ans) = summary$Experimental.Method
round(ans,2)
```

    ##               X-Ray                 NMR Electron Microscopy               Other 
    ##               89.06                8.13                2.51                0.19 
    ##        Multi Method 
    ##                0.10

What propotions are protein?

``` r
round(sum(summary$Proteins) / sum(summary$Total),2)
```

    ## [1] 0.93

``` r
library(bio3d)
pdb = read.pdb("1hsg")
```

    ##   Note: Accessing on-line PDB file

``` r
pdb
```

    ## 
    ##  Call:  read.pdb(file = "1hsg")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

``` r
aa321(pdb$atom[,"resid"])
```

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH
    
    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ##    [1] "P" "P" "P" "P" "P" "P" "P" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "I" "I"
    ##   [19] "I" "I" "I" "I" "I" "I" "T" "T" "T" "T" "T" "T" "T" "L" "L" "L" "L" "L"
    ##   [37] "L" "L" "L" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "Q"
    ##   [55] "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "R" "R" "R" "R" "R" "R" "R" "R" "R" "R"
    ##   [73] "R" "P" "P" "P" "P" "P" "P" "P" "L" "L" "L" "L" "L" "L" "L" "L" "V" "V"
    ##   [91] "V" "V" "V" "V" "V" "T" "T" "T" "T" "T" "T" "T" "I" "I" "I" "I" "I" "I"
    ##  [109] "I" "I" "K" "K" "K" "K" "K" "K" "K" "K" "K" "I" "I" "I" "I" "I" "I" "I"
    ##  [127] "I" "G" "G" "G" "G" "G" "G" "G" "G" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q"
    ##  [145] "L" "L" "L" "L" "L" "L" "L" "L" "K" "K" "K" "K" "K" "K" "K" "K" "K" "E"
    ##  [163] "E" "E" "E" "E" "E" "E" "E" "E" "A" "A" "A" "A" "A" "L" "L" "L" "L" "L"
    ##  [181] "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "D" "D" "D" "D" "D" "D" "D"
    ##  [199] "D" "T" "T" "T" "T" "T" "T" "T" "G" "G" "G" "G" "A" "A" "A" "A" "A" "D"
    ##  [217] "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "T" "T" "T"
    ##  [235] "T" "T" "T" "T" "V" "V" "V" "V" "V" "V" "V" "L" "L" "L" "L" "L" "L" "L"
    ##  [253] "L" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E"
    ##  [271] "E" "M" "M" "M" "M" "M" "M" "M" "M" "S" "S" "S" "S" "S" "S" "L" "L" "L"
    ##  [289] "L" "L" "L" "L" "L" "P" "P" "P" "P" "P" "P" "P" "G" "G" "G" "G" "R" "R"
    ##  [307] "R" "R" "R" "R" "R" "R" "R" "R" "R" "W" "W" "W" "W" "W" "W" "W" "W" "W"
    ##  [325] "W" "W" "W" "W" "W" "K" "K" "K" "K" "K" "K" "K" "K" "K" "P" "P" "P" "P"
    ##  [343] "P" "P" "P" "K" "K" "K" "K" "K" "K" "K" "K" "K" "M" "M" "M" "M" "M" "M"
    ##  [361] "M" "M" "I" "I" "I" "I" "I" "I" "I" "I" "G" "G" "G" "G" "G" "G" "G" "G"
    ##  [379] "I" "I" "I" "I" "I" "I" "I" "I" "G" "G" "G" "G" "G" "G" "G" "G" "F" "F"
    ##  [397] "F" "F" "F" "F" "F" "F" "F" "F" "F" "I" "I" "I" "I" "I" "I" "I" "I" "K"
    ##  [415] "K" "K" "K" "K" "K" "K" "K" "K" "V" "V" "V" "V" "V" "V" "V" "R" "R" "R"
    ##  [433] "R" "R" "R" "R" "R" "R" "R" "R" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Y"
    ##  [451] "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "D" "D" "D" "D" "D" "D" "D"
    ##  [469] "D" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "I" "I" "I" "I" "I" "I" "I" "I"
    ##  [487] "L" "L" "L" "L" "L" "L" "L" "L" "I" "I" "I" "I" "I" "I" "I" "I" "E" "E"
    ##  [505] "E" "E" "E" "E" "E" "E" "E" "I" "I" "I" "I" "I" "I" "I" "I" "C" "C" "C"
    ##  [523] "C" "C" "C" "G" "G" "G" "G" "H" "H" "H" "H" "H" "H" "H" "H" "H" "H" "K"
    ##  [541] "K" "K" "K" "K" "K" "K" "K" "K" "A" "A" "A" "A" "A" "I" "I" "I" "I" "I"
    ##  [559] "I" "I" "I" "G" "G" "G" "G" "T" "T" "T" "T" "T" "T" "T" "V" "V" "V" "V"
    ##  [577] "V" "V" "V" "L" "L" "L" "L" "L" "L" "L" "L" "V" "V" "V" "V" "V" "V" "V"
    ##  [595] "G" "G" "G" "G" "P" "P" "P" "P" "P" "P" "P" "T" "T" "T" "T" "T" "T" "T"
    ##  [613] "P" "P" "P" "P" "P" "P" "P" "V" "V" "V" "V" "V" "V" "V" "N" "N" "N" "N"
    ##  [631] "N" "N" "N" "N" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I"
    ##  [649] "I" "I" "G" "G" "G" "G" "R" "R" "R" "R" "R" "R" "R" "R" "R" "R" "R" "N"
    ##  [667] "N" "N" "N" "N" "N" "N" "N" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L"
    ##  [685] "L" "L" "L" "L" "L" "T" "T" "T" "T" "T" "T" "T" "Q" "Q" "Q" "Q" "Q" "Q"
    ##  [703] "Q" "Q" "Q" "I" "I" "I" "I" "I" "I" "I" "I" "G" "G" "G" "G" "C" "C" "C"
    ##  [721] "C" "C" "C" "T" "T" "T" "T" "T" "T" "T" "L" "L" "L" "L" "L" "L" "L" "L"
    ##  [739] "N" "N" "N" "N" "N" "N" "N" "N" "F" "F" "F" "F" "F" "F" "F" "F" "F" "F"
    ##  [757] "F" "P" "P" "P" "P" "P" "P" "P" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "I"
    ##  [775] "I" "I" "I" "I" "I" "I" "I" "T" "T" "T" "T" "T" "T" "T" "L" "L" "L" "L"
    ##  [793] "L" "L" "L" "L" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W" "W"
    ##  [811] "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "R" "R" "R" "R" "R" "R" "R" "R" "R"
    ##  [829] "R" "R" "P" "P" "P" "P" "P" "P" "P" "L" "L" "L" "L" "L" "L" "L" "L" "V"
    ##  [847] "V" "V" "V" "V" "V" "V" "T" "T" "T" "T" "T" "T" "T" "I" "I" "I" "I" "I"
    ##  [865] "I" "I" "I" "K" "K" "K" "K" "K" "K" "K" "K" "K" "I" "I" "I" "I" "I" "I"
    ##  [883] "I" "I" "G" "G" "G" "G" "G" "G" "G" "G" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q"
    ##  [901] "Q" "L" "L" "L" "L" "L" "L" "L" "L" "K" "K" "K" "K" "K" "K" "K" "K" "K"
    ##  [919] "E" "E" "E" "E" "E" "E" "E" "E" "E" "A" "A" "A" "A" "A" "L" "L" "L" "L"
    ##  [937] "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L" "D" "D" "D" "D" "D" "D"
    ##  [955] "D" "D" "T" "T" "T" "T" "T" "T" "T" "G" "G" "G" "G" "A" "A" "A" "A" "A"
    ##  [973] "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "D" "T" "T"
    ##  [991] "T" "T" "T" "T" "T" "V" "V" "V" "V" "V" "V" "V" "L" "L" "L" "L" "L" "L"
    ## [1009] "L" "L" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E" "E"
    ## [1027] "E" "E" "M" "M" "M" "M" "M" "M" "M" "M" "S" "S" "S" "S" "S" "S" "L" "L"
    ## [1045] "L" "L" "L" "L" "L" "L" "P" "P" "P" "P" "P" "P" "P" "G" "G" "G" "G" "R"
    ## [1063] "R" "R" "R" "R" "R" "R" "R" "R" "R" "R" "W" "W" "W" "W" "W" "W" "W" "W"
    ## [1081] "W" "W" "W" "W" "W" "W" "K" "K" "K" "K" "K" "K" "K" "K" "K" "P" "P" "P"
    ## [1099] "P" "P" "P" "P" "K" "K" "K" "K" "K" "K" "K" "K" "K" "M" "M" "M" "M" "M"
    ## [1117] "M" "M" "M" "I" "I" "I" "I" "I" "I" "I" "I" "G" "G" "G" "G" "G" "G" "G"
    ## [1135] "G" "I" "I" "I" "I" "I" "I" "I" "I" "G" "G" "G" "G" "G" "G" "G" "G" "F"
    ## [1153] "F" "F" "F" "F" "F" "F" "F" "F" "F" "F" "I" "I" "I" "I" "I" "I" "I" "I"
    ## [1171] "K" "K" "K" "K" "K" "K" "K" "K" "K" "V" "V" "V" "V" "V" "V" "V" "R" "R"
    ## [1189] "R" "R" "R" "R" "R" "R" "R" "R" "R" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q"
    ## [1207] "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y" "D" "D" "D" "D" "D" "D"
    ## [1225] "D" "D" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "Q" "I" "I" "I" "I" "I" "I" "I"
    ## [1243] "I" "L" "L" "L" "L" "L" "L" "L" "L" "I" "I" "I" "I" "I" "I" "I" "I" "E"
    ## [1261] "E" "E" "E" "E" "E" "E" "E" "E" "I" "I" "I" "I" "I" "I" "I" "I" "C" "C"
    ## [1279] "C" "C" "C" "C" "G" "G" "G" "G" "H" "H" "H" "H" "H" "H" "H" "H" "H" "H"
    ## [1297] "K" "K" "K" "K" "K" "K" "K" "K" "K" "A" "A" "A" "A" "A" "I" "I" "I" "I"
    ## [1315] "I" "I" "I" "I" "G" "G" "G" "G" "T" "T" "T" "T" "T" "T" "T" "V" "V" "V"
    ## [1333] "V" "V" "V" "V" "L" "L" "L" "L" "L" "L" "L" "L" "V" "V" "V" "V" "V" "V"
    ## [1351] "V" "G" "G" "G" "G" "P" "P" "P" "P" "P" "P" "P" "T" "T" "T" "T" "T" "T"
    ## [1369] "T" "P" "P" "P" "P" "P" "P" "P" "V" "V" "V" "V" "V" "V" "V" "N" "N" "N"
    ## [1387] "N" "N" "N" "N" "N" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I" "I"
    ## [1405] "I" "I" "I" "G" "G" "G" "G" "R" "R" "R" "R" "R" "R" "R" "R" "R" "R" "R"
    ## [1423] "N" "N" "N" "N" "N" "N" "N" "N" "L" "L" "L" "L" "L" "L" "L" "L" "L" "L"
    ## [1441] "L" "L" "L" "L" "L" "L" "T" "T" "T" "T" "T" "T" "T" "Q" "Q" "Q" "Q" "Q"
    ## [1459] "Q" "Q" "Q" "Q" "I" "I" "I" "I" "I" "I" "I" "I" "G" "G" "G" "G" "C" "C"
    ## [1477] "C" "C" "C" "C" "T" "T" "T" "T" "T" "T" "T" "L" "L" "L" "L" "L" "L" "L"
    ## [1495] "L" "N" "N" "N" "N" "N" "N" "N" "N" "F" "F" "F" "F" "F" "F" "F" "F" "F"
    ## [1513] "F" "F" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1531] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1549] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1567] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1585] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1603] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1621] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1639] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1657] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"
    ## [1675] "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X" "X"

``` r
protein_select = atom.select(pdb,"protein", value = TRUE)
write.pdb(protein_select, "protein_select.pdb")
```

``` r
ligand_select = atom.select(pdb, "ligand", value = TRUE)
write.pdb(ligand_select, "ligand_select.pdb")
ligand_select
```

    ## 
    ##  Call:  trim.pdb(pdb = pdb, sele)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 45,  XYZs#: 135  Chains#: 1  (values: B)
    ## 
    ##      Protein Atoms#: 0  (residues/Calpha atoms#: 0)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 45  (residues: 1)
    ##      Non-protein/nucleic resid values: [ MK1 (1) ]
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

``` r
#library(bio3d.view)
#view(ligand_select)
```
