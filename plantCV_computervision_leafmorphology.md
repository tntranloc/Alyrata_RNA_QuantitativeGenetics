




   


```python
import plantcv as pcv
import cv2
import numpy as np
import os
```


```python
from plantcv import plantcv as pcv
from plantcv.parallel import WorkflowInputs
```


```python
# Input/output options
args = WorkflowInputs(
    images=["A.jpg"],    
    names="R",
    result="R_res",
    outdir=".",
    writeimg=True,
    debug="plot",
    sample_label="genotype"
    )
```


```python
# Set debug to the global parameter 
pcv.params.debug = args.debug

# Set plotting size (default = 100)
pcv.params.dpi = 100

# Increase text size and thickness to make labels clearer
# (size may need to be altered based on original image size)
pcv.params.text_size = 10
pcv.params.text_thickness = 20
```


```python
# Read in your image, which is based on the path you put above
img, path, filename = pcv.readimage(filename="/Users/nhutran/Documents/PhD/epidom/assays/leaf_photo2/R.jpg")
```


    
![png](output_11_0.png)
    



```python
# Crop image if necessary. This is optional. 
crop_img = pcv.crop(img=img, x=0, y=0, h=3500, w=4500)
```


    
![png](output_12_0.png)
    



```python
# We define the circular region interest by x,y for the center, and r for the radius of the circle to get made 
reference_roi = pcv.roi.circle(img=crop_img, x=150, y=450, r=100)
```


    
![png](output_13_0.png)
    



```python
scaled_img = cv2.resize(crop_img, None, fx=1, fy=1)
```


```python
# no need grey scaled as I scanned the leaf herbarium and the image is also very close to grey scale
# if not, make sure to convert it
# gray_img = pcv.rgb2gray(scaled_img)
```


    
![png](output_15_0.png)
    



```python
# no need binary conversion if image is already grey scale
#binary_threshold_img = pcv.threshold.binary(gray_img, threshold=120)

```


    
![png](output_16_0.png)
    



```python
#segmented_img, obj = pcv.morphology.segment_skeleton(skel_img=binary_threshold_img)

```


    
![png](output_17_0.png)
    



```python
#leaf_obj, other_obj = pcv.morphology.segment_sort(skel_img=binary_threshold_img,
                                                  objects=obj)
```


    
![png](output_18_0.png)
    



```python
#colorspaces = pcv.visualize.colorspaces(rgb_img=scaled_img, original_img= True)

```


```python
# take channel L here as it visualise the leaves the best
a = pcv.rgb2gray_lab(rgb_img=scaled_img, channel='l')
```


    
![png](output_20_0.png)
    



```python
# adjust the threshold to see it best, higher threshold means less background noise/sharper image
a_thresh = pcv.threshold.binary(gray_img=a, threshold=150, object_type='dark')

```


![png](output_21_0.png)
    

```python
# adjust the size here similarly
a_fill = pcv.fill(bin_img=a_thresh, size=120)
```
    
![png](output_22_0.png)
    

```python
# define areas of interest
rois = pcv.roi.auto_grid(mask=a_fill, nrows=4, ncols=6, img=crop_img)

```
    
![png](output_23_0.png)
    

```python
labeled_mask, num_plants = pcv.create_labels(mask=a_fill, rois=rois, roi_type="partial")
 
```


    
![png](output_24_0.png)
    



```python
shape_img = pcv.analyze.size(img=crop_img, labeled_mask=labeled_mask, n_labels=24)
 
```

   
![png](output_25_0.png)
    


```python
pcv.outputs.save_results(filename="/Users/nhutran/Documents/PhD/epidom/assays/leaf_photo2/result/R")
```

