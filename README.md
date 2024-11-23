Materials for the paper: Kowshika Sarker, Ruoqing Zhu, Hannah D. Holscher, and ChengXiang Zhai. 2023. Augmenting nutritional metabolomics with a genome-scale metabolic model for assessment of diet intake. In Proceedings of the 14th ACM International Conference on Bioinformatics, Computational Biology, and Health Informatics (BCB '23). Association for Computing Machinery, New York, NY, USA, Article 4, 1–10. https://doi.org/10.1145/3584371.3612958
- ```preprocess-metabolome.py``` Filters raw baseline and end metabolite concentrations to remove samples measured at only one timestamp, metabolites with too many missing values, and imputes the remaining missing values.
- ```preprocess-human-gem.py``` Filters Human-GEM to produce the 9 reaction sets.
- ```compute-change-feature.py``` Calculates the proposed 'Change' features with preprocessed metabolome and a given reaction set.
- ```compute-ratio-feature.py``` Calculates the proposed 'Ratio' features with preprocessed metabolome and a given reaction set.
- ```compute-prob-feature.py``` Calculates the proposed 'Prob' features with preprocessed metabolome and a given reaction set.
- ```classification.py``` Performs hyperparameter tuning and classification using random forests with a given feature matrix.
- ```preprocess-and-compute-features.py``` Preprocesses metabolomic profiles and Human-GEM, computes features integrating these two processed sources.
