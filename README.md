# Using prior knowledge for channel compression in highly multiplexed imaging

This repository contains the code used to obtain the results presented in the thesis. It contains two directories, one for each analyzed dataset.

## IMC dataset
The code for preprocessing the IMC images and the complete workflow for channel compression are in `workflow.R`. The statistics and figures presented in the results section can be generated using the `results.R` script. The raw images can be found [here](https://data.mendeley.com/datasets/cydmwsfztj/2).

## CODEX dataset
The preprocessed images can be downloaded by running the `download.sh` script. For segmentation, [CellSeg](https://github.com/michaellee1/CellSeg) was called with the parameters specified in `batch_segmentation.py`. Running `rename_channels.sh` reorganizes the segmented images and standardizes the channel names between batches. The `merge_data.R` script aggregates the data into a SpatialExperiment object. The entire channel compression pipeline is contained in `workflow.R`. The statistics and figures were generated with `results.R`. To run [STELLAR](https://github.com/snap-stanford/stellar), two of the files from the default installation must be replaced. The modified files are in the `stellar` directory. They change the input and output format, not the algorithm itself.