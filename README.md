 # BriLo
This repository provides tools to download PSP/WISPR Level 3 data, select and track coronal rays, fit their measured brightness over elongation with a theoretical profile of the Thomson Scattering theory to infer positional and physical properties of coronal streamer observed by PSP.

First clone the repository:

    cd <directory_to_save_repository>
    git clone https://github.com/gretacappello/BriLo
    cd BriLo
   
Then install requirements using the environment.yml file:

    conda env create -f environment.yml
    conda activate brilo


Tracking
-----
Go to the directory in which the codes are stored and open the first code which is a Jupyter Notebook:

    jupyter notebook downloadWISPR_filter_point_track_save.ipynb

Follow the notebook steps to track the selected ray and store the coordinates along the ray.

Notes:
- Make sure you select ray and pointy rays. Avoid blurry and over structured rays.
- You may want to repeat the track at least three times to ensure to have a good extimation of the errors.

Fitting
-----
Once the ray has been tracked, you can use the brilo.py script to find the optimal solution for C and γ based on the selected measurement points. To run the Python script, execute:

     python brilo.py

Before running, check and edit the initial settings in brilo.py to ensure you have specified the correct example name, the path to the L3 WISPR file, the path to the directory where figures should be saved and any other relevant keywords or configuration options. Make sure these settings match your tracked ray and data files to obtain accurate results.


Contact and referencing
-----
A study on coronal rays observed by WISPR using BriLo is currently under development.

For any question or request please contact me at: greta.cappello@uni-graz.at 

Please cite {...} paper and the code when used for science pubblication.

The version v1.0 of BriLo is archived on Zenodo at DOI {...}


