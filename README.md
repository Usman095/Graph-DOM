# Graph-DOM
Graph-DOM is a tool to analyze comprehensive fragmentation data for dissolved organic matter (DOM) using graph algorithms.  

To setup and use Graph-DOM, follow the instructions below:

## Setup
1. Ubuntu 16.04 or later.
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
10. Set parameters in `config.ini` file. See the Config section for details about setting parameters.
