# Force-Inference for epithelial cells

`@Author: Ali Hashmi`

see my Wolfram Community post for more information: https://community.wolfram.com/groups/-/m/t/1571507

This is a Mathematica implementation of the `Bayesian Force Inference` technique by Ishihara and Sugimura (Journal of Theoretical Biology, 2012). The script file `Force Inference.m` will infer tension between epithelial cells and the pressure within the cells by using a binarized image of an epithelia.

##### Please note that the borders of the image used below for the demo, `image.tif` (inside for ReadMe folder), have pixels that are filled. 

https://www.biorxiv.org/content/early/2018/11/23/475012 (for more details, check recent work by Lenne Lab)

###### The notebook (.nb file) contains the same code as the script file. The only difference is that I have included some auxiliary parts in the notebook to check the robustness of detection of vertices and some visualizations to ensure the implementation is fool-proof. The extra results are included in `robustness check.pdf`



![alt-text](https://github.com/alihashmiii/Force-Inference/blob/master/for%20ReadMe/im1.png)


![alt-text](https://github.com/alihashmiii/Force-Inference/blob/master/for%20ReadMe/im2.png)


![alt-text](https://github.com/alihashmiii/Force-Inference/blob/master/for%20ReadMe/im3.png)


![alt-text](https://github.com/alihashmiii/Force-Inference/blob/master/for%20ReadMe/im4.png)



#### Robustness Check
<a href="https://github.com/alihashmiii/Force-Inference/blob/master/robustness%20check.pdf"><img src= "https://github.com/alihashmiii/Force-Inference/blob/master/for%20ReadMe/robustnesscheckImg.png" alt="Illustration" width="400px"/></a> &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp;
