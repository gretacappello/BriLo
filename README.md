# BriLo

This repository provides tools to download WISPR Level 3 (L3) data, select and track coronal rays, and fit their measured brightness over elongation using a theoretical profile based on Thomson Scattering theory. This allows you to infer the positional and physical properties of coronal streamers observed by PSP/WISPR.

## Getting Started

1. **Clone the repository**
   ```bash
   cd <directory_to_save_repository>
   git clone https://github.com/gretacappello/BriLo
   cd BriLo

    Install requirements

conda env create -f environment.yml
conda activate wispr-ray-tracker

Open the Jupyter Notebook

    jupyter notebook downloadWISPR_filter_point_track_save.ipynb

    Follow the notebook steps to:

        Download WISPR Level 3 files.

        Select a coronal ray of interest.

        Track points along the ray.

        Fit the tracked points with a quadratic function.

        Save the fitted coordinates for later analysis.

Notes

    Make sure you select well-defined, pointed rays. Avoid blurry or highly structured rays.

    It is recommended to repeat the tracking process at least three times to obtain a reliable estimate of the errors.

Running the Fitting Script

Once the ray has been tracked, you can use the brilo.py script to find the optimal solution for CC and γγ based on the selected measurement points. To run the script, execute:

python brilo.py

Before running, check and edit the initial settings in brilo.py to ensure you have specified:

    The correct example name.

    The path to the L3 WISPR file.

    The path to the directory where figures should be saved.

    Any other relevant keywords or configuration options.

Make sure these settings match your tracked ray and data files to obtain accurate results.
