# Graph-DOM
Graph-DOM is a tool to analyze comprehensive fragmentation data for dissolved organic matter (DOM) using graph algorithms.  

To setup and use Graph-DOM, follow the instructions below:

## Setup
1. System with Ubuntu 16.04 or later with at least 8 CPU cores and 120GBs of memory is recommended.
2. Install Anaconda using these [instructions](https://docs.conda.io/en/latest/miniconda.html#linux-installers).
3. [Download](https://github.com/Usman095/Graph-DOM/archive/refs/tags/v1.0.tar.gz) the source code.
4. Extract the `tar.gz` file:  
   `tar -xzf Graph-DOM-1.0.tar.gz`
6. Change the current working directory to the downloaded Graph-DOM directory in the previous step. 
7. Create virtual environment with project dependencies using the following command:  
   `conda create --name graph-dom --file requirements.txt python=3.9`
8. Type `y` and press enter when prompted `Proceed ([y]/n)?`.
9. Activate the virtual environment:
   `conda activate graph-dom`
10. Set parameters in `config.ini` file. See the [Config](#config) section for details about setting parameters.
11. Run Graph-DOM using the following command: `python3 main.py`
12. The program will generate pathways for each precursor and then analyze the pathways to generate families. Once completed, the output file and plots can located in the `output` directory. For details on the output files and plots, see the publication `Place holder for publication.`


## Config
`config.ini` file contains user adjustable parameters to run Graph-DOM.
1. `num_cores`: Number of cores to be used for generating pathways in parallel. Should be <= number of CPU cores of the system being used.
2. `use_NS`: Whether to use Nitrogen and Sulphur as elements for creating pathways.
3. `input_file_path`: Relative path to the input file. Input files are provided inside the Input directory.
4. `multiple`: Multiple of each neutral loss to consider to find the next peak.
5. `tolerance` Fragment tolerance for generating pathways.
6. `nominal_tolerance`: Fragments within +-nominal_tolerance Da will be considered precursors.
7. `overlap_len`: Overlap length threshold of pathways for creating families. Should be at least 2. To enforce complete overlap set it to a large number e.g. 100
