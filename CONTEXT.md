# ECODA Domain Glossary

## Benchmark cohort
A biological scRNA-seq dataset included in the cross-cohort embedding method comparison. The verification-only debug subset is not a publication cohort.

## Cell type annotation method
A method that assigns cell-level biological labels or operational cluster identities that can define cell-composition categories. Author labels, automated labels, and Leiden clusters are annotation sources.

## Cell-type identity
A valid biological label or intentional operational cluster category assigned by an annotation source. Missing, blank, and explicit placeholder values such as `NA`, `None`, or `Unknown` are not identities.

## Cell-type count
The number of distinct valid cell-type identities available to one annotation method within one benchmark cohort. For inferred composition methods, the corresponding inferred component categories are counted instead.

## Leiden cluster
An unsupervised cluster identity used as an operational annotation category. A Leiden cluster is not automatically a biological cell type, even when it is displayed under the cell-type count metric.

## Deconvolution component
An inferred cell-composition category produced without direct per-cell annotation. Its categories are method outputs, not source cell-level identities.

## Embedding method
A method that represents samples in a feature or distance space.
