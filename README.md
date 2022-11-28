# AGREDA_1.1

New repository of genome-scale metabolic models of the human gut microbiota, developed starting from AGREDA repository (https://github.com/tblasco/AGREDA). This new version of the metabolic network add 80 new phenolic compounds and their related 195 reactions predicted using an enzyme promiscuity algorithm called RetroPath RL (Koch, M., Duigou, T. & Faulon, J. L. Reinforcement learning for bioretrosynthesis. ACS Synth. Biol. 9, 157–168 (2020)).

AGREDA_1.1 mixed-bag model and subsequent species models can be found at AGREDA_v1.1.0.zip file.

For further information, please refer to:

Francesco Balzerani fbalzerani@tecnun.es
Xabier Cendoya xcendoya@tecnun.es
Telmo Blasco tblasco@tecnun.es
Iñigo Apaolaza iaemparanza@tecnun.es
Francisco J. Planes fplanes@tecnun.es

Citing AGREDA_1.1

Balzerani, Francesco, et al. "Prediction of degradation pathways of phenolic compounds in the human gut microbiota through enzyme promiscuity methods." NPJ systems biology and applications 8.1 (2022): 1-9.

PIPELINE
Please refer to the following sections for a pipeline execution description.

SOFTWARE REQUIREMENTS
It will be necessary to install Matlab. We encourage users to download release R2018a.

Free academic licenses for the IBM CPLEX solver can be obtained from https://www.ibm.com/developerworks/community/blogs/jfp/entry/CPLEX_Is_Free_For_Students?lang=en. We encourage users to download CPLEX version 12.8.0.

Cobratoolbox github link https://github.com/opencobra/cobratoolbox.

It will be necessary to install Python. We encourage users to download version >=3.6.

EXECUTION
The execution is divided in two steps, as shown below. Within each step, various scripts can be run indipendetly since in the input and output folder are present all the necessary files.

01_FilesGenerationAndRetroPathRLApplication: Files generation and application of the enzyme promiscuity algorithm RetroPath RL.
02_ManualRevisionAndReconstruction: Extraction of the results produced by the algorithm, their analysis and revision followed by the metabolic network reconstruction

OUTPUT

In 02_ManualRevisionAndReconstruction/output/outputReconstruction there are two type of outputs:
  - AGREDA1_1.mat: supra-organism model
  - singleSpecies folder: contains the reconstruction for the 818 organisms considering the improvements made with AGREDA1_1 (they are compatible with cobrapy).



