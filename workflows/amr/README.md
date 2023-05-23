## CZ ID AMR Workflow

CZ ID's AMR pipeline implements the [Resistance Gene Identifier (RGI)](https://github.com/arpcard/rgi) tool for AMR sequence detection. The RGI tool is used to compare quality controlled reads and assembled contigs against AMR references sequences from the [Comprehensive Antibiotic Resistance Database](https://card.mcmaster.ca/) (CARD).  Further documentation on how to use the CZ ID AMR workflow, including a [pipeline workflow](https://chanzuckerberg.zendesk.com/hc/en-us/articles/15091031482644-AMR-Pipeline-Workflow) can be found in the [CZ ID help center](https://chanzuckerberg.zendesk.com/hc/en-us/categories/15001531592980-Antimicrobial-Resistance-Analysis).

## Changelog

v1.2.3 -- May 24, 2023 -- Initial pipeline release

## Reference Files

Filename | Provenance 
---------|-----------
s3://czid-public-references/card/2023-05-22/* | All files downloaded from https://card.mcmaster.ca/download on May 23, 2023. CARD version: 3.2.6 Wildcard version: 4.0.0
