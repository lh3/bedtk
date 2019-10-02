Bedtk is a set of simple tools to process BED files. It so far implements
intersection, sorting, merging and computing the breadth of coverage. Bedtk is
not as versatile as [bedtools][bedtools] and never aims to match the bedtools
feature set. It is instead developed mainly for performance focused tasks. It
is several to tens of times faster and uses a fraction of memory. Bedtk also
provides a few convenient functions. For example, sorting, merging and
intersection can be done in one go without Unix pipes.

[bedtools]: https://github.com/arq5x/bedtools2
[cr]: https://github.com/lh3/cgranges
