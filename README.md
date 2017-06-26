# Vectorized Stress Majorization
Source code of paper "Revisiting Stress Majorization as a Unified Framework for
Interactive Constrained Graph Visualization" by Yunhai Wang, Yanyan Wang et al.

Project page: ...

## Note
Due to the complexity of our user interactive system, the source code of the interaction part is not included in our distribution. The UI system is provided as an executable file for Windows operating system.

These codes include the main contribution of our algorithm, Vectorized Stress Majorization(VSM). We provide two sample codes for you to test. You can modify them to fit your need.

Note that the edges, distances matrix, circles, and communities are read from files. We do not provide the algorithms for calculating distances or finding communities and circles. And the format of data should be organized as the sample datasets we provide.

We modified viennaCL(http://viennacl.sourceforge.net/) to enhance the performance for our need particularly.

## Compilation

We provided CMake lists for both *nix and Win32 environments and Microsoft Visual Studio projects for Windows.

Linux:
Use CMakeLists_NIX.txt for CMake.
CUDA SDK required.

Windows:
Use CMakeLists.txt with CMake 3.8 or later for MSVC on windows.
Visual Studio 2015 solution file is also provided.

Configurations are tested on Ubuntu 16.04 with Clang and Visual Studio 2015 on Windows.

## Dependencies
Our system is running under CUDA 8.0, viennaCL, Eigen, and Glut in visual studio 2015. If you want to change the configuration of our solution file, you can modify .... and reload the solution file in visual studio.

## Cite
You can use our codes for research purpose only. And please cite our paper when you use our codes.
@article{wang2018vsm,
  title={ Revisiting Stress Majorization as a Unified Framework for Interactive Constrained Graph Visualization},
  author={ Yunhai Wang, Yanyan Wang, Yinqi Sun, Lifeng Zhu, Kecheng Lu, Chi-Wing Fu, Michael Sedlmair, Oliver Deussen, and Baoquan Chen,
  journal={IEEE Trans. Vis. & Comp. Graphics (Proc. IEEE Information Visualization (Infovis) 2017)},
  year={2017},
  publisher={IEEE}
}


## Licence
Vectorized Stress Majorization is open-sourced software licensed under the MIT license.

## Contact
If you find any bugs or have any ideas of optimizing these codes, please contact me via yanyanwang93 [at] gmail [dot] com.

Yanyan Wang
