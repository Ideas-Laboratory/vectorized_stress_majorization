# Vectorized Stress Majorization
Source code of paper "Revisiting Stress Majorization as a Unified Framework for Interactive Constrained Graph Visualization" by Yunhai Wang and Yanyan Wang et al. conditionally accepted for publication in IEEE Information Visualization (InfoVis) 2017.

Project page: ...

## Note
Due to the complexity of our user interactive system, the source code of the interaction part is not included in our distribution, but we provide an executable file of the UI system: run on 64-bit Windows 10 operating system. 

The distribution includes our algorithm, the Vectorized Stress Majorization (VSM). In addition, we provide two sample code for testing. You can modify the code based on your needs.

Note that the edges, distances matrix, circles, and communities are read from files. We do not provide the algorithms for calculating distances or finding communities and circles.  The format of the data should be organized as the sample datasets we provided.

Lastly, we modified viennaCL (http://viennacl.sourceforge.net/) to improve the performance for our need particularly. The modified viennaCL is provided in vectorized_stress_majorization/src/VectorizedStressMajorization/viennacl.

## Compilation
We provided CMake lists for both *nix and Win32 environments and Microsoft Visual Studio projects for Windows.

•	Linux: Use CMakeLists_NIX.txt for CMake. CUDA SDK required.

•	Windows: Use CMakeLists.txt with CMake 3.8 or later for MSVC on windows. Visual Studio 2015 solution file is also provided.
Configurations are tested on Ubuntu 16.04 with Clang and Visual Studio 2015 on Windows.

## Dependencies
Our system runs on CUDA 8.0, viennaCL, Eigen, and GLUT in Visual Studio 2015. 

## Reference
Our code is provided for research purpose only.
Please cite our paper when you use our codes.

@article{wang2018vsm,

  title={Revisiting Stress Majorization as a Unified Framework for Interactive Constrained Graph Visualization},
  
  author={Yunhai Wang, Yanyan Wang, Yinqi Sun, Lifeng Zhu, Kecheng Lu, Chi-Wing Fu, Michael Sedlmair, Oliver Deussen, and Baoquan Chen},
  
  journal={IEEE Transactions Visualization and Computer Graphics (IEEE Information Visualization (InfoVis) 2017)},
  
  year={2017},
  
  publisher={IEEE}
  
}


## Licence
Vectorized Stress Majorization is an open-sourced software licensed under the MIT license.

## Contact
If you find any bugs or have any ideas of optimizing these codes, please contact me via cloudseawang [at] gmail [dot] com.

Yanyan Wang
