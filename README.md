# BriLo
This repository provides tools to download WISPR Level 3 (L3) data, select and track coronal rays, fit their measured brightness over elongation with a theoretical profile of the Thomson Scattering theory to infer positional and physical properties of coronal streamer observed by PSP/WISPR

1. Clone the repository:
    cd <directory_to_save_repository>
    git clone https://github.com/gretacappello/BriLo
    cd BriLo
   
3. Install requirements:

    conda env create -f environment.yml
   
5. Open the Jupyter Notebook:

    jupyter notebook downloadWISPR_filter_point_track_save.ipynb

7. Follow the notebook steps to:

- Download WISPR Level 3 files.
- Select a coronal ray of interest.
- Track points along the ray.
- Fit the tracked points with a quadratic function.
- Save the fitted coordinates for later analysis.

Notes
-----
- Make sure you select ray and pointy rays. Avoid blurry and over structured rays.
- You may want to repeat the track at least three time to ensure to have a good extimation of the errors.

Once the ray has been tracked, you can use the brilo.py script to find the optimal solution for C and γ based on the selected measurement points. To run the script, execute:

python brilo.py

Before running, check and edit the initial settings in brilo.py to ensure you have specified: 
The correct example name.
The path to the L3 WISPR file.
The path to the directory where figures should be saved.
Any other relevant keywords or configuration options.

Make sure these settings match your tracked ray and data files to obtain accurate results.

  
