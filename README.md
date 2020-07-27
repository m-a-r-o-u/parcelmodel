# A atmospheric parcel model

This model was developed to investigate the impact of thermal radiation and saturation fluctuations due to turbulence
on the development of the cloud droplet distribution. The motivation behind it, is to improve our understanding
of rain development in ice free clouds. The description and results are summarized in the research paper:

[Broadening of the Cloud Droplet Size Distribution due to Thermal Radiative Cooling: Turbulent Parcel Simulations](https://journals.ametsoc.org/jas/article/77/6/1993/345221)

## Basic Usage

Install the [pipenv](https://pipenv-fork.readthedocs.io/en/latest/install.html) package manager:
```
pip install --user pipenv
```
Clone the repository and install the required python packages:
```
pipenv install
```
Startup the pipenv shell
```
pipenv shell
```
Create the config.yaml input file in the root folder:
```
initial_conditions:
    w: 1
    l: 50
    feedback:
        latent_heat: 1
        radiation: 1
    dt: 1
    t_max: 1200
    dt_output: 20
    T: 300
    p: 100000
    qv: 17.e-3
    z0: 1000
    particle_distribution:
        type: double_log_normal
        ratio: 0.6
        mu1: 1.e-8
        sigma1: 1
        mu2: 1.e-7
        sigma2: 1
        total: 1.e+8
        Nsd: 100
    radiation_schema:
        type: stefan_boltzmann
        l: 50
        factor: 1
    turbulence_schema:
        type: markov_schema
        Nsd: 100
        l: 50
        epsilon: 50.e-6
        dt: 1
    atmosphere_schema:
        type: linear_and_hydrostatic
output:
    logger:
        type: NetCDFLogger

globals:
    executer: numpy
    file_path: ./output/
```
And finally, run the following command to create the results as an netcdf file in ./output/:
```
./run.py config.yaml
```
