```python
import plantcv as pcv 
import pandas as pd
import os
import numpy as np
import cv2
```


```python
from plantcv import plantcv as pcv
from plantcv.parallel import WorkflowInputs
```


```python
import os
```


```python
######### My Batch 1 ########
```


```python
os.chdir("/Users/nhutran/Documents/batch1")
```


```python
### batch 1- different time points ###
```


```python
#os.chdir("/Users/nhutran/Documents/batch1/24.10")
#os.chdir("/Users/nhutran/Documents/batch1/30.10")
#os.chdir("/Users/nhutran/Documents/batch1/06.11")
os.chdir("/Users/nhutran/Documents/batch1/13.11")
```


```python
# Input/output options
args = WorkflowInputs(
    images=["G.JPG"],    
    names="G_1311",
    result="G_1311_res",
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
# Inputs:
#   filename = Image file to be read in 
#   mode     = How to read in the image; either 'native' (default), 
#              'rgb', 'gray', 'csv', or 'envi'
img, path, filename = pcv.readimage(filename="/Users/nhutran/Documents/batch1/13.11/G.JPG")
```


    
![png](output_9_0.png)
    



```python
# We define the circular region interest by x,y for the center, and r for the radius of the circle to get made 
reference_roi = pcv.roi.circle(img=img, x=1375, y=1360, r=370)
```


    
![png](output_10_0.png)
    



```python
#real-life reference object's size # here is diameter of the pots
pot_dia_mm = 90
# Calculate the scaling factor based on the diameter of the pots
pixels_per_mm = (370/2)/(pot_dia_mm/2) # 180 is radius of the reference_roi = r/2 = 360/2
scaling_factor = 1/pixels_per_mm  # Invert the ratio to get the scaling factor

```


```python
# Resize the image based on the scaling factor
scaled_img = cv2.resize(img, None, fx=scaling_factor, fy=scaling_factor)
```


```python
# Continue with PlantCV analysis on the scaled image

```


```python
colorspaces = pcv.visualize.colorspaces(rgb_img=scaled_img, original_img=False)

```


    
![png](output_14_0.png)
    



```python
a = pcv.rgb2gray_lab(rgb_img=scaled_img, channel='a')
```


    
![png](output_15_0.png)
    



```python
a_thresh = pcv.threshold.binary(gray_img=a, threshold=110, object_type='dark')

```


    
![png](output_16_0.png)
    



```python
a_fill = pcv.fill(bin_img=a_thresh, size=10)
```


    
![png](output_17_0.png)
    



```python
rois = pcv.roi.auto_grid(mask=a_fill, nrows=4, ncols=6, img=scaled_img)

```


    
![png](output_18_0.png)
    



```python
# Create a labeled mask, this function works very similarly to the roi.filter step above 
labeled_mask, num_plants = pcv.create_labels(mask=a_fill, rois=rois, roi_type="partial")

```


    
![png](output_19_0.png)
    



```python
shape_img = pcv.analyze.size(img=scaled_img, labeled_mask=labeled_mask, n_labels=24)

```


    
![png](output_20_0.png)
    



```python

```


```python
pcv.outputs.save_results(filename='/Users/nhutran/Documents/batch1/b1_results/G_1311')
```


```python

```


```python
######### My Batch 2 ########
```


```python
os.chdir("/Users/nhutran/Documents/batch2")
```


```python
### batch 2- different time points ###
```


```python
os.chdir("/Users/nhutran/Documents/batch2/25.10")
#os.chdir("/Users/nhutran/Documents/batch2/31.10")
#os.chdir("/Users/nhutran/Documents/batch2/06.11")
#os.chdir("/Users/nhutran/Documents/batch2/13.11")
#os.chdir("/Users/nhutran/Documents/batch2/20.11")
```


```python
# Input/output options
args = WorkflowInputs(
    images=["L.JPG"],    
    names="L_2510",
    result="L_2510_res",
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
# Inputs:
#   filename = Image file to be read in 
#   mode     = How to read in the image; either 'native' (default), 
#              'rgb', 'gray', 'csv', or 'envi'
img, path, filename = pcv.readimage(filename="/Users/nhutran/Documents/batch2/25.10/L.JPG")
```


    
![png](output_30_0.png)
    



```python
# We define the circular region interest by x,y for the center, and r for the radius of the circle to get made 
reference_roi = pcv.roi.circle(img=img, x=1630, y=1770, r=450)
```


    
![png](output_31_0.png)
    



```python
#real-life reference object's size # here is diameter of the pots
pot_dia_mm = 90
# Calculate the scaling factor based on the diameter of the pots
pixels_per_mm = (450/2)/(pot_dia_mm/2) # 180 is radius of the reference_roi = r/2 = 360/2
scaling_factor = 1/pixels_per_mm  # Invert the ratio to get the scaling factor

```


```python
# Resize the image based on the scaling factor
scaled_img = cv2.resize(img, None, fx=scaling_factor, fy=scaling_factor)
```


```python
colorspaces = pcv.visualize.colorspaces(rgb_img=scaled_img, original_img=False)

```


    
![png](output_34_0.png)
    



```python
a = pcv.rgb2gray_lab(rgb_img=scaled_img, channel='a')
```


    
![png](output_35_0.png)
    



```python
a_thresh = pcv.threshold.binary(gray_img=a, threshold=107, object_type='dark')

```


    
![png](output_36_0.png)
    



```python
a_fill = pcv.fill(bin_img=a_thresh, size=10)
```


    
![png](output_37_0.png)
    



```python
rois = pcv.roi.auto_grid(mask=a_fill, nrows=3, ncols=5, img=scaled_img)

```


    
![png](output_38_0.png)
    



```python
# Create a labeled mask, this function works very similarly to the roi.filter step above 
labeled_mask, num_plants = pcv.create_labels(mask=a_fill, rois=rois, roi_type="partial")

```


    
![png](output_39_0.png)
    



```python
shape_img = pcv.analyze.size(img=scaled_img, labeled_mask=labeled_mask, n_labels=15)

```


    
![png](output_40_0.png)
    



```python
pcv.outputs.save_results(filename='/Users/nhutran/Documents/batch2/b2_results/L_2510')
```


```python

```


```python
#### My batch 3 - different time points ####
```


```python
#os.chdir("/Users/nhutran/Documents/batch3/06.11")
#os.chdir("/Users/nhutran/Documents/batch3/13.11")
os.chdir("/Users/nhutran/Documents/batch3/20.11")
```


```python
# Input/output options
args = WorkflowInputs(
    images=["O.JPG"],    
    names="O_2011",
    result="O_2011_res",
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
# Inputs:
#   filename = Image file to be read in 
#   mode     = How to read in the image; either 'native' (default), 
#              'rgb', 'gray', 'csv', or 'envi'
img, path, filename = pcv.readimage(filename="/Users/nhutran/Documents/batch3/20.11/O.JPG")
```


    
![png](output_47_0.png)
    



```python
# We define the circular region interest by x,y for the center, and r for the radius of the circle to get made 
reference_roi = pcv.roi.circle(img=img, x=1480, y=1350, r=400)
```


    
![png](output_48_0.png)
    



```python
#real-life reference object's size # here is diameter of the pots
pot_dia_mm = 90
# Calculate the scaling factor based on the diameter of the pots
pixels_per_mm = (400/2)/(pot_dia_mm/2) # 180 is radius of the reference_roi = r/2 = 360/2
scaling_factor = 1/pixels_per_mm  # Invert the ratio to get the scaling factor

```


```python
# Resize the image based on the scaling factor
scaled_img = cv2.resize(img, None, fx=scaling_factor, fy=scaling_factor)
```


```python
colorspaces = pcv.visualize.colorspaces(rgb_img=scaled_img, original_img=False)

```


    
![png](output_51_0.png)
    



```python
a = pcv.rgb2gray_lab(rgb_img=scaled_img, channel='a')
```


    
![png](output_52_0.png)
    



```python
a_thresh = pcv.threshold.binary(gray_img=a, threshold=111, object_type='dark')

```


    
![png](output_53_0.png)
    



```python
a_fill = pcv.fill(bin_img=a_thresh, size=10)
```


    
![png](output_54_0.png)
    



```python
rois = pcv.roi.auto_grid(mask=a_fill, nrows=3, ncols=6, img=scaled_img)

```


    
![png](output_55_0.png)
    



```python
labeled_mask, num_plants = pcv.create_labels(mask=a_fill, rois=rois, roi_type="partial")

```


    
![png](output_56_0.png)
    



```python
shape_img = pcv.analyze.size(img=scaled_img, labeled_mask=labeled_mask, n_labels=24)

```


    
![png](output_57_0.png)
    



```python
pcv.outputs.save_results(filename='/Users/nhutran/Documents/batch3/b3_results/O_2011')
```


```python

```


```python

```


```python

```
