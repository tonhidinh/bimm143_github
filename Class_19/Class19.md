# Class 19: Mini Project Cancer Mutation
Tonhi Dinh (A18508020)

I downloaded my sequence from the class website and moved it to my
working directory:

From a **BLASTp** searach against human refseq I see that it is:

> Q1. \[1pt\] What protein do these sequences correspond to? (Give both
> full gene/protein name and official symbol).

- Official symbol: BRAF
- Official name: B-Raf proto-oncogene, serine/threonine kinase

``` r
library(bio3d)

seq <- read.fasta("./Class 19/A18508020_mutant_seq.fa")
seq
```

                   1        .         .         .         .         .         60 
    wt_healthy     MAALSGGGGGGAEPGQALFNGDMEPEAGAGAGAAASSAADPAIPEEVWNIKQMIKLTQEH
    mutant_tumor   MAALSGGGGGGAEPGQALFNGDMEPEAGAGAGAAASSAADPAIPEEVWNIKQMIKLTQEH
                   ************************************************************ 
                   1        .         .         .         .         .         60 

                  61        .         .         .         .         .         120 
    wt_healthy     IEALLDKFGGEHNPPSIYLEAYEEYTSKLDALQQREQQLLESLGNGTDFSVSSSASMDTV
    mutant_tumor   IEALLDKFGGEHNPPSIYLEAYEEYTSKLDALQQREQQLLESLGNGTDFSVSSSASMDTV
                   ************************************************************ 
                  61        .         .         .         .         .         120 

                 121        .         .         .         .         .         180 
    wt_healthy     TSSSSSSLSVLPSSLSVFQNPTDVARSNPKSPQKPIVRVFLPNKQRTVVPARCGVTVRDS
    mutant_tumor   TSSSSSSLSVLPSSLSVFQNPTDVARSNPKSPQKPIVRVFLPNKQRTVVPARCGVTVRDS
                   ************************************************************ 
                 121        .         .         .         .         .         180 

                 181        .         .         .         .         .         240 
    wt_healthy     LKKALMMRGLIPECCAVYRIQDGEKKPIGWDTDISWLTGEELHVEVLENVPLTTHNFVRK
    mutant_tumor   LKKALMMRGLIPECCAVYRIQDGEKKPIGWDTDISWLTGEELHVEVLENVPLTTHNFVRK
                   ************************************************************ 
                 181        .         .         .         .         .         240 

                 241        .         .         .         .         .         300 
    wt_healthy     TFFTLAFCDFCRKLLFQGFRCQTCGYKFHQRCSTEVPLMCVNYDQLDLLFVSKFFEHHPI
    mutant_tumor   TFFTLAFCDFCRKLLFQGFRCQTCGYKFHQRCSTEVPLMCVNYDQLDLLFVSKFFEHHPI
                   ************************************************************ 
                 241        .         .         .         .         .         300 

                 301        .         .         .         .         .         360 
    wt_healthy     PQEEASLAETALTSGSSPSAPASDSIGPQILTSPSPSKSIPIPQPFRPADEDHRNQFGQR
    mutant_tumor   PQEEASLAETALTSGSSPSAPASDSIGPQILTSPSPSKSIPIPQPFRPADEDHRNQFGQR
                   ************************************************************ 
                 301        .         .         .         .         .         360 

                 361        .         .         .         .         .         420 
    wt_healthy     DRSSSAPNVHINTIEPVNIDDLIRDQGFRGDGGSTTGLSATPPASLPGSLTNVKALQKSP
    mutant_tumor   DRSSSAPNVHINTIEPVNIDDLIRDQGFRGDGGSTTGLSATPPASLPGSLTNVKALQKSP
                   ************************************************************ 
                 361        .         .         .         .         .         420 

                 421        .         .         .         .         .         480 
    wt_healthy     GPQRERKSSSSSEDRNRMKTLGRRDSSDDWEIPDGQITVGQRIGSGSFGTVYKGKWHGDV
    mutant_tumor   GPQRERKSSSSSEDRNRMKTLGRRDSSDDWEIPDGQITVGQRIGSGSFGTVYKGKWHGDV
                   ************************************************************ 
                 421        .         .         .         .         .         480 

                 481        .         .         .         .         .         540 
    wt_healthy     AVKMLNVTAPTPQQLQAFKNEVGVLRKTRHVNILLFMGYSTKPQLAIVTQWCEGSSLYHH
    mutant_tumor   AVKMLNVTAPTPQQLQAFKNEVGVLRKTRHVNILLFMGYSTKPQLAIVTQWCEGSSLYHH
                   ************************************************************ 
                 481        .         .         .         .         .         540 

                 541        .         .         .         .         .         600 
    wt_healthy     LHIIETKFEMIKLIDIARQTAQGMDYLHAKSIIHRDLKSNNIFLHEDLTVKIGDFGLATV
    mutant_tumor   LHIIETKFEMIKLIDIARQTAQGMDYLHAKSIIHRDLVSNNIFLHEDLTVEIGDFGLATV
                   ************************************* ************ ********* 
                 541        .         .         .         .         .         600 

                 601        .         .         .         .         .         660 
    wt_healthy     KSRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAFGIVLYELMTGQLPYSNIN
    mutant_tumor   KSRWSGSHQFEQLSGSILWRAPEVIRMQDKNPYSFQSDVYAFGIVLYELMTGQLPYSNIN
                   ******************* **************************************** 
                 601        .         .         .         .         .         660 

                 661        .         .         .         .         .         720 
    wt_healthy     NRDQIIFMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELLARS
    mutant_tumor   NRDQIYFMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELLARS
                   ***** ****************************************************** 
                 661        .         .         .         .         .         720 

                 721        .         .         .         .     766 
    wt_healthy     LPKIHRSASEPSLNRAGFQTEDFSLYACASPKTPIQAGGYGAFPVH
    mutant_tumor   LPKIHRSASEPSLNRAGFQTEDFSLYACASPKTPIQAGGYGAFPVH
                   ********************************************** 
                 721        .         .         .         .     766 

    Call:
      read.fasta(file = "./Class 19/A18508020_mutant_seq.fa")

    Class:
      fasta

    Alignment dimensions:
      2 sequence rows; 766 position columns (766 non-gap, 0 gap) 

    + attr: id, ali, call

We could score residue conservation and then find the non 1.0 scoring
position. These will be the mutation positions:

``` r
conserv(seq)
```

      [1]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
     [16]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
     [31]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
     [46]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
     [61]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
     [76]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
     [91]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [106]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [121]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [136]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [151]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [166]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [181]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [196]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [211]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [226]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [241]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [256]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [271]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [286]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [301]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [316]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [331]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [346]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [361]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [376]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [391]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [406]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [421]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [436]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [451]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [466]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [481]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [496]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [511]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [526]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [541]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [556]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [571]  1.0  1.0  1.0  1.0  1.0  1.0  1.0 -0.3  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [586]  1.0  1.0  1.0  1.0  1.0  0.1  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [601]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [616]  1.0  1.0  1.0  1.0 -0.2  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [631]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [646]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [661]  1.0  1.0  1.0  1.0  1.0 -0.2  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [676]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [691]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [706]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [721]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [736]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [751]  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
    [766]  1.0

``` r
mutation.sites <- which(conserv(seq) < 1)
```

> Q2. \[6pts\] What are the tumor specific mutations in this particular
> case ( e.g. A130V)?

``` r
paste(seq$ali[1, mutation.sites],
mutation.sites,
seq$ali[2, mutation.sites])
```

    [1] "K 578 V" "K 591 E" "M 620 R" "I 666 Y"

> Q3. \[1pts\] Do your mutations cluster to any particular domain and if
> so give the name and PFAM id of this domain? Alternately note whether
> your protein is single domain and provide it’s PFAM id/accession and
> name (e.g. PF00613 and PI3Ka).

Domains from hmmer search

![](Screenshot%202026-03-10%20at%2011.44.23.png)

PFAM:

Description:

Coordinates: - K578V: Transferase(Phosphotransferase) domain 1 - K591E:
Transferase(Phosphotransferase) domain 1 - M620R:
Transferase(Phosphotransferase) domain 1 - I666Y:
Transferase(Phosphotransferase) domain 1

> Q4. \[2pts\] Using the NCI-GDC list the observed top 2 missense
> mutations in this protein (amino acid substitutions)?

V640E - chr7:g.140753336A\>T V640M - chr7:g.140753337C\>T

> Q5. \[2pts\] What two TCGA projects have the most cases affected by
> mutations of this gene? (Give the TCGA “code” and “Project Name” for
> example “TCGA-BRCA” and “Breast Invasive Carcinoma”).

TCGA-THCA - Thyroid Carcinoma TCGA-SKCM - Skin Cutaneous Melanoma

> Q6. \[3pts\] List one RCSB PDB identifier with 100% identity to the
> wt_healthy sequence and detail the percent coverage of your query
> sequence for this known structure? Alternately, provide the most
> similar in sequence PDB structure along with it’s percent identity,
> coverage and E-value. Does this structure “cover” (i.e. include or
> span the amino acid residue positions) of your previously identified
> tumor specific mutations?

Chain B, Serine/threonine-protein kinase B-raf \[Homo sapiens\]

- Query coverage: 39%

- Percent identity: 100.00%

- E-value: 0.0

![](4MNE%20copy.png)

> Q7. \[10pts\] Using AlphaFold notebook generate a structural model
> using the default parameters for your mutant sequence.

> Q8. \[2pts\] Considering only your mutations in high quality structure
> regions (with a pLDDT score \> 70) are any of the mutations on the
> surface of the protein and hence have a potential to interfere with
> protein-protein interaction events? List these mutations below
> (e.g. A130V)

> Q9. \[5pts\] Please comment on how useful and/or reliable you think
> your AlphaFold structural model is for your entire sequence and the
> main domain where your mutations lie? You may wish to compare your
> model to the PDB structure you found in Q6.

> Q10. \[10pts\] Are any of the identified “hot spots” near your cancer
> specific mutation sites or the most commonly mutated sites from the
> NCI-GDC? If so which mutation site(s)?
