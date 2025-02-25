## Code examples for hu.MAP 3.0 training, testing, and generating predictions

All code is available as Jupyter notebooks

### Accessing feature matrix files

All feature matrices can be pulled from our [dataset](https://huggingface.co/datasets/sfisch/hu.MAP3.0) on HuggingFace
and examples of using the *datasets* module can be seen in **hu.MAP3.0_testing.ipynb** and **hu.MAP3.0_training.ipynb**

In brief:
  ```python
  from autogluon.tabular import TabularPredictor
  from datasets import load_dataset

  dataset = load_dataset('sfisch/hu.MAP3.0')
  train = dataset["train"].to_pandas()
  test = dataset["test"].to_pandas()
  ```

The full feature matrix can be pulled using the *huggingface_hub* module as seen in **generating_predictions_w_hu.MAP3.0.ipynb** or accessed via the downloads page of [our website](https://humap3.proteincomplexes.org)

To pull from HuggingFace and use as a pandas dataframe:
  ```python
  from huggingface_hub import hf_hub_download
  import pandas as pd

  full_file = hf_hub_download(repo_id="sfisch/hu.MAP3.0", filename='full/humap3_full_feature_matrix_20220625.csv.gz', repo_type='dataset')
  full_featmat = pd.read_csv(full_file, compression="gzip")
  ```
### Accessing the hu.MAP 3.0 model

The [hu.MAP 3.0 model](https://huggingface.co/sfisch/hu.MAP3.0_AutoGluon) can also be downloaded from HuggingFace using the *huggingface_hub* module as seen in **generating_predictions_w_hu.MAP3.0.ipynb** and **hu.MAP3.0_testing.ipynb**. Note, that version 0.4.0 of AutoGluon is needed to use the model.

  ```python
  from autogluon.tabular import TabularPredictor
  from huggingface_hub import snapshot_download

  model_dir = snapshot_download(repo_id="sfisch/hu.MAP3.0_AutoGluon")
  predictor = TabularPredictor.load(f"{model_dir}/huMAP3_20230503_complexportal_subset10kNEG_notScaled_accuracy")
  ```
### Accessing other testing/training data

All test/train complexes and interactions can be found on [our website](https://humap3.proteincomplexes.org/download) or on [HuggingFace](https://huggingface.co/datasets/sfisch/hu.MAP3.0/tree/main/reference_interactions)
