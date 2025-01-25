This repository contains the code used for the publication of Tsiros et al. (2025) entitled " Optimizing Nanoparticle-Mediated Drug Delivery: Insights from Compartmental Modeling via the CompSafeNano Cloud Platform".

* The "Nanoparticle_delivery_application_publication.R" contains the code used for generating the results presented in the publication.

* File "model_fitting.R" can be used to fit the model to new data. Specifically, based on the examples.csv file, the user is requested to provide NP concentration-time
profiles for the NP under study for the corresponding model compartments, i.e., the administration site, cell vicinity, cell interior, off-target sites and excreta.
The time vector should be provided in hours and the masses in μg. Note that all files should be in the same working directory, which also needs to be set as the working
R directory. Finally, end-users can predefine some parameter values at specific values, according to the description of the code. Once set, the whole script should be
excecuted. Upon completion, the optimal parameter values and goodness-of-fit plots will be generated

* Having run the "model_fitting.R" file, end-users can then run the "local_sensitivity.R" code to obtain the local sensitivity of NP delivery to the cell interior with
respect to time.
