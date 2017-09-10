This work is released under GPLv3 license as in GPLv3.pdf and GPLv3.tex with copyright as in copyright.txt

For proper functionality, please add all folders and subfolders of this base_pc_v1 folder to MATLAB path.

basis has functions related to basis identification and evaluation
examples has some ready made examples to help understand how some functions are used
license has GPLv3 license data as also found in this folder
other has functions that are not readily categorized into the other folders
qoi_evals has different functions qoi evaluations
sampler has functions associated with the random sampling, coherence-optimal sampling, and sample adaptivity
solver has functions related to surrogate identification

default_parameters.m in other is a function that may be useful to tailor to the machine on which this code is running, or to the specific problem being investigated.
Similarly, example_5.m and example_6.m are coded with parameters that are defined without the use of default_parameters.m.

To cite this software, please cite the following paper.

@article{base_pc_2017,
  author    = {J. Hampton and A. Doostan},
  journal = {ArXiv e-prints},
  title     = {Basis Adaptive Sample Efficient Polynomial Chaos (BASE-PC)},
  year      = {2017},
  note      = {Available from https://arxiv.org/abs/1702.01185}
}
 
To cite coherence-optimal sampling, whose functionality is included, please cite the following papers.
@article{CohOptSamplingEll1,
  title={Compressive Sampling of Polynomial Chaos Expansions: Convergence Analysis and Sampling Strategies},
  author={J. Hampton and A. Doostan},
  journal={Journal of Computational Physics},
  volume={280},
  pages={363--386},
  year={2015}
}


@article{CohOptSamplingEll2,
  title={Coherence motivated sampling and convergence analysis of least squares polynomial Chaos regression},
  author={J. Hampton and A. Doostan},
  journal={Computer Methods in Applied Mechanics and Engineering},
  volume={290},
  pages={73--97},
  year={2015},
  publisher={Elsevier}
}