# mutrate-mlp
Feed-forward neural net to predict cancer mutation rate from histone modifications

# reproducibility 

* run all cells in `data_prep.ipynb`  
  * this will download and prepare data from the Epigenetic Roadmap Project, and The Cancer Genome Atlas
  * should complete in under 2 min with good network connection

* in `train.ipynb`, modify `cancer` and `feature_subset` as desired

* run all cells in `train.ipynb`
  * this will train the selected model, and prepare plots of test performance, feature interpretation
  * completion time depends on your compute
