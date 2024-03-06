# Augmenting nutritional metabolomics with a genome-scale metabolic model for assessment of diet intake
## Scripts
| Filename  | Functionality |
| ------------- | ------------- |
|preprocess-metabolome.py| Filters raw baseline and end metabolite concentrations to remove samples measured at only one timestamp, metabolites with too many missing values, and imputes the remaining missing values.|
|preprocess-human-gem.py| Filters Human-GEM to produce the 9 reaction sets.|
|compute-change-feature.py|Calculates the proposed 'Change' features with preprocessed metabolome and a given reaction set.|
|compute-ratio-feature.py|Calculates the proposed 'Ratio' features with preprocessed metabolome and a given reaction set.|
|compute-prob-feature.py|Calculates the proposed 'Prob' features with preprocessed metabolome and a given reaction set.|
|classification.py|Performs hyperparameter tuning and classification using random forests with a given feature matrix.|
| preprocess-and-compute-features.py | Preprocesses metablomic profiles and Human-GEM, computes features integerating these two processed sources.|
