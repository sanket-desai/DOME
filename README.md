# DOmain Mutation Estimator (DOME) - user manual

DOME, developed by [Dutt lab](http://www.actrec.gov.in/pi-webpages/AmitDutt/dutt_index.html), [ACTREC-TMC](https://actrec.gov.in/))

An algorithm to prioritize low-frequency somatic mutations integrating somatic mutation propensity, mutation distribution across domain and functional and biochemical residue context within protein structures.

## Getting Started

DOME is developed using python3 and R shiny package. The database is provided as a Tabix indexed (.tbi) TSV file along with the toolkit package for download.

### Pre-requisites required for installation of DOME

Users are required to download DOME toolkit package from the [webpage](http://www.actrec.gov.in/pi-webpages/AmitDutt/DOME/DOME.html). Approximately 2 GB of system space is required to store the database toolkit package. The downloaded package needs to be extracted using "tar" or "winrar" extractor, depending on the operating system used. The package contains the following contents:

```terminal
dome (root folder)
    ├── Manual (User Manual)
    ├── app.R (Rshiny App)
    ├── data
    ├── input_dome.txt (example: input file)
    ├── output
    ├── domeenv.yml
    └── src (DOME algorithm and accesory scripts)

```

**System Prerequisites:**

  - [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)

The installation of the required packages is performed through Anaconda/Miniconda, which needs to be installed. With Conda installed, users can use the following commands in the Linux / Unix environment to install the required python packages.

### Installation

DOME can be downloaded from this [webpage](http://www.actrec.gov.in/pi-webpages/AmitDutt/DOME/DOME.html). Untar the downloaded tar.gz file using the following command:
```bash
tar -xvzf dome_v0.2.tar.gz
cd dome
```
Make sure that the complete 'untarred' data directory is present in the dome home directory before proceeding with the further installation and running dome.
In the home directory of dome (dome), run the following command:

```bash
conda env create -f domeenv.yml
```
This installs all the dependencies and tools required to run DOME.

For running any script pertaining to DOME, make sure that the domeenv Conda environment is activated. More information about Conda environments can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

### Running DOME toolkit

The required scripts to run DOME toolkit are placed in the “src” directory.

**To load the R Shiny based GUI or via web-browser use the following command: **

**In Linux terminal:**

```bash
$ conda activate domeenv
$ R -e "shiny::runApp('app.R', launch.browser = TRUE)"
```

![](/Screenshots/Execution.png)

**In Windows CMD**:

```bash
> "Path/to/R.exe" -e "shiny::runApp("app.R", launch.browser = TRUE)"
```



### Guide to use DOME toolkit

Once users executes the above commands DOME toolkit GUI appears.

Following are the detailed description for each of the tabs-

1. File Upload

   ![](/Screenshots/Final_Dome_Score.png)

   

### DOME output

Gene based or mutation (position) based search returns a table containing the following columns:

| Column  names                 | Description                                                  |
| ----------------------------- | ------------------------------------------------------------ |
| Entry name                    | Entry name of protein  as per Uniprot database               |
| Position                      | Amino acid position  in the protein                          |
| Reference                     | Reference amino acid  in the protein                         |
| Altered                       | Altered amino acid in  the protein                           |
| Domain                        | Domain name                                                  |
| Mutation type                 | H - hotspot, R -  resistant, '-' - not present in significant mutation database |
| Analogous to                  | Query protein  position is analogous to this known significant mutation |
| CADD                          | CADD score (converted  functional score, as described in manuscript) |
| COSMIC count                  | Mutation count in the  COSMIC database                       |
| ClinVar                       | Mutation reported in  ClinVar                                |
| Site                          | A functional site as  described in Uniprot database          |
| Phosphosite                   | A site defined in  PhosphoSite database                      |
| Proximity to  significant     | Geometrical (3D)  proximity to a significant mutation (hotspot/resistant) |
| Proximity to  functional site | Geometrical (3D)  proximity to a functional site described in Uniprot |
| Proximity to  PSP site        | Geometrical (3D)  proximity to a site defined in PhosphoSite database |
| Proximity to  3D hotspot      | Geometrical (3D)  proximity to a defined 3D hotspot (hotspotdb) |
| Is interface                  | Residue is a part of  a defined protein interface (3D PIP database) |
| Entropy                       | Entropy score (refer  methodology) of the residue within domain |
| DOME Score                    | Score computed  using the DOME algorithm scoring metric (scaled 0-1) |
