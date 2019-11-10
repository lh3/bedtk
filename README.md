## Getting Started

```sh
# compile and install
git clone https://github.com/lh3/bedtk
cd bedtk && make
# filter a BED or VCF file
./bedtk flt test/test-anno.bed.gz test/test-iso.bed.gz
./bedtk flt -cw100 test/test-anno.bed.gz test/test-sub.vcf.gz
# intersect (no sorting needed; overlapping records allowed)
./bedtk isec test/test-anno.bed.gz test/test-iso.bed.gz
# compute the breadth of coverage
./bedtk cov test/test-anno.bed.gz test/test-iso.bed.gz
# sort a BED file
./bedtk sort test/test-iso.bed.gz
./bedtk sort -s test/chr_list.txt test/test-iso.bed.gz
# merge overlapping records (no sorting needed)
./bedtk merge test/test-anno.bed.gz
```

## Introduction

Bedtk is a set of simple tools to process BED files. It so far implements
intersection, subtraction, sorting, merging and computing the breadth of
coverage. Bedtk is not as versatile as [bedtools][bedtools] and never aims to
match the bedtools feature set. It instead focuses on performance. Bedtk
is several to tens of times faster and uses a fraction of memory. It also
provides a few convenient functions. For example, sorting, merging and
intersection can be done in one go without Unix pipes.

[bedtools]: https://github.com/arq5x/bedtools2
[cr]: https://github.com/lh3/cgranges
