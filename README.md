# Dicrimination-and-learning-of-temporal-input-sequences-in-a-cerebellar-Purkinje-cell-model

These programs were used in the article "Discrimination and learning of temporal input sequences in a cerebellar Purkinje cell model". 
(https://www.frontiersin.org/articles/10.3389/fncel.2023.1075005/full)

In this research, I modified part of the cerebellar circuit model (Kobayashi et al. 2021), so I only provide the modified program. Please download the original model by the author's github (https://github.com/TairaKobayashi/Cerebellum-neuron-model-simulation-using-RKC) and replace the files. Note that the program "stdp.c" in chapter3-3 is independent from the cerebellar model and was used to generate stimulation values after learning.

For running these programs, you need three input files ( input_time.txt, input_pos.txt, iput_size.txt ). This Repository have sample files.

Kobayashi, T., Kuriyama, R. & Yamazaki, T. Testing an Explicit Method for Multi-compartment Neuron Model Simulation on a GPU. Cogn Comput (2021). https://doi.org/10.1007/s12559-021-09942-6
