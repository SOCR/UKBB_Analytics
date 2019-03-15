# UKBB_Analytics
**SOCR UK Biobank Data Analytics (Mental health and mood disorders)**

This [SOCR](http://socr.umich.edu/) [GitHUb](https://github.com/SOCR) partition includes the end-to-end computational protocol, results, 
validation and scripts supporting a Big Data discovery study involving the [UK Biobank data](http://www.ukbiobank.ac.uk).


# Overview
The [UK Biobank data](http://www.ukbiobank.ac.uk) represents an extremely rich health research resource. It provides enormous
opportunities for data scientists, computational analysts and bioinformaticians to examine, model, and analyze census-like multisource healthcare data. This archive presents many difficult challenges related to aggregation and harmonization of complex data,
feature heterogeneity, and robust health analytics. 

In this study, we used 7,614 imaging, clinical, and phenotypic features of 9,914 subjects to perform deep computed phenotyping using machine learning and artificial intelligence techniques. We determined the top 20 most salient features contributing to the separation of all participants into separate clusters (cohorts). By jointly representing and modeling the significant clinical and demographic variables along with the derived salient neuroimaging features we used decision rules to predict the presence and progression of depression. 

This study reported on the consistency and reliability of the derived computed phenotypes and the top salient imaging biomarkers that
contributed to the unsupervised clustering. The main outcome is a clinical decision support system that identifies and jointly utilizes the most critical biomarkers for predicting mental health, e.g., depression. The study externally validated the method by applying the technique to an independent (out-of-bag) testing dataset. Such approaches may lead to reducing healthcare expenses and improving the processes of diagnosis, forecasting, and tracking of normal and pathological aging.


# Resources
* All **data** used in this study is openly and publicly available from the [UK Biobank data](http://www.ukbiobank.ac.uk). Access requires appropriate.
* The **code directory** contains the complete collection of R scripts used to preprocess, model, visualize, interrogate and understand the data.
* The **results** folder includes some of the derived measures and result summaries.


# Usage
In the spirit of [open science](https://en.wikipedia.org/wiki/Open_science) and in accordance with the [SOCR Licensing](http://socr.umich.edu/html/SOCR_CitingLicense.html), all materials here are [CC-BY](https://creativecommons.org/licenses/) and [LGPL](https://opensource.org/licenses/lgpl-license) licensed.

# Contact
[SOCR Personnel](http://www.socr.umich.edu/people/).


# Reference
* Zhou, Y, Zhao, l, Zhou, N, Zhao, Y, Marino, S, Wang, T, Sun, H, Toga, AW, [Dinov, ID](http://www.socr.umich.edu/people/dinov/).  (2019) *Predictive Big Data Analytics using the UK Biobank Data*, 
[Scientific Reports](https://www.nature.com/srep/), *in press*, [DOI : 10.1038/s41598-019-41634-y](http://doi.org/10.1038/s41598-019-41634-y).
