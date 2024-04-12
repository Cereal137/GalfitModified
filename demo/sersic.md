# Sérsic Index

At local universe, Sérsic Index increase with observed wavelength because bulges (large n) are more dominant at longer wavelength. Also, differential dust attenuation can account for larger Sérsic Index at longer wavelength.

However, in my GALFITM output, there seems no tendency toward increasing Sérsic index with wavelength. In my GALFITM configuration, I make position angle, axis ratio of light profile to be congruent between all bands, and leave all band magnitude and centroid position to be varying independently. The centroid position are required to be only slightly deviating from the initial guess. The Sérsic index and effective radius are limited by a linear function to ensure a monotonic trend.

![](slope_hist.png)

## Overall Trend

The overall trend is shown below. Their rest-frame wavelengths are infered from EAZY photometric redshift.

![](sersic_rf.png)
![](sersic_rf2.png)

The growth of Sérsic Index towards redder end is hard to perceive.

The counter-intuitive situation could be either the case that our GALFITM configuration failed to represent the truth, or that at higher redshift, stellar population of bulges could be less redder than expected (intrinsically or reddened by concentrated dust), which requires further tests.

## Detailed Analysis

Notaion:
Label 'Linear' denotes the GALFITM configuration, where the Sérsic Index and effective radius can only be varying linearly with wavelength. Label 'Single Band' represents GALFIT results from single band images, where there's no constraint between bands. And results from linear configuration are shown in the left column, while results from single band images are shown in the right column.

By visually inspecting the images, I divided the sample into 2 main catagories: isolated samples and blended samples. Blended samples are sources surrounded by one or more bright companions, which require multi-component fitting. For those samples, I apply additional Sérsic components to the model.

![](magdist.png)

## Examples of Increasing Sérsic Index

Galaxies with increasing Sérsic Index are slightly smaller in number when compared to decreasing ones. Their distribution in magnitude is shown below.

![](magdist_inc.png)

Samples that are defined as isolated are dominant in those which have increasing Sérsic index. However, the difference cannot account for the neutral trend in Sérsic Index, which is illustrated in the catagorized trend shown above.

Here are some examples of increasing Sérsic Index:

![](sersicplot/31079.png)
<div>
<img src="galfitmimage_linear/nircam1_31079.png" alt="" style="width: 49%; display: inline-block;">
<img src="singlebandimage/nircam1_31079.png" alt="" style="width: 49%; display: inline-block;">
</div>

![](sersicplot/24329.png)
<div>
<img src="galfitmimage_linear/nircam1_24329.png" alt="" style="width: 49%; display: inline-block;">
<img src="singlebandimage/nircam1_24329.png" alt="" style="width: 49%; display: inline-block;">
</div>

![](sersicplot/38682.png)
<div>
<img src="galfitmimage_linear/nircam1_38682.png" alt="" style="width: 49%; display: inline-block;">
<img src="singlebandimage/nircam1_38682.png" alt="" style="width: 49%; display: inline-block;">
</div>

![](sersicplot/40785.png)
<div>
<img src="galfitmimage_linear/nircam1_40785.png" alt="" style="width: 49%; display: inline-block;">
<img src="singlebandimage/nircam1_40785.png" alt="" style="width: 49%; display: inline-block;">
</div>

![](sersicplot/37221.png)
<div>
<img src="../multicomp/img_mc_multi/37221.png" alt="" style="width: 49%; display: inline-block;">
<img src="../multicomp/img_mc_single/nircam6_37221.png" alt="" style="width: 49%; display: inline-block;">
</div>


<br/><br/><br/><br/><br/><br/><br/><br/><br/><br/>

## Examples of decreasing Sérsic Index

It's more important to take a serious look at samples with decreasing Sérsic Index. Their distribution in magnitude is shown below.

![](magdist_dec.png)

More blended samples are found in those with decreasing Sérsic index, but not very significantly.

<br/><br/><br/><br/><br/><br/><br/><br/><br/><br/><br/><br/><br/><br/>
<br/><br/><br/><br/><br/><br/><br/><br/><br/>

Here are some examples that could have decreasing Sérsic Index. The results are categorized by their mean magnitude, which serves as proxy for their image quality, or SNR. Samples with brighter mean magnitude are expected to be robust in single band Sérsic Index measurement (GALFIT), if 1D profile fitting is reliable.

### 1. Mean Mag < 25 mag

![](sersicplot/23225.png)
<div>
<img src="galfitmimage_linear/nircam1_23225.png" alt="" style="width: 49%; display: inline-block;">
<img src="singlebandimage/nircam1_23225.png" alt="" style="width: 49%; display: inline-block;">
</div>

![](sersicplot/25075.png)
<div>
<img src="galfitmimage_linear/nircam1_25075.png" alt="" style="width: 49%; display: inline-block;">
<img src="singlebandimage/nircam1_25075.png" alt="" style="width: 49%; display: inline-block;">
</div>

![](sersicplot/12831.png)
<div>
<img src="galfitmimage_linear/nircam6_12831.png" alt="" style="width: 49%; display: inline-block;">
<img src="singlebandimage/nircam6_12831.png" alt="" style="width: 49%; display: inline-block;">
</div>

![](sersicplot/18581.png)
<div>
<img src="../multicomp/img_mc_multi/18581.png" alt="" style="width: 49%; display: inline-block;">
<img src="../multicomp/img_mc_single/nircam6_18581.png" alt="" style="width: 49%; display: inline-block;">
</div>

### 2. Mean Mag > 25 mag but < 26.5 mag

![](sersicplot/26553.png)
<div>
<img src="galfitmimage_linear/nircam1_26553.png" alt="" style="width: 49%; display: inline-block;">
<img src="singlebandimage/nircam1_26553.png" alt="" style="width: 49%; display: inline-block;">
</div>

![](sersicplot/30649.png)
<div>
<img src="galfitmimage_linear/nircam1_30649.png" alt="" style="width: 49%; display: inline-block;">
<img src="singlebandimage/nircam1_30649.png" alt="" style="width: 49%; display: inline-block;">
</div>

![](sersicplot/40073.png)
<div>
<img src="galfitmimage_linear/nircam1_40073.png" alt="" style="width: 49%; display: inline-block;">
<img src="singlebandimage/nircam1_40073.png" alt="" style="width: 49%; display: inline-block;">
</div>

![](sersicplot/20851.png)
<div>
<img src="galfitmimage_linear/nircam2_20851.png" alt="" style="width: 49%; display: inline-block;">
<img src="singlebandimage/nircam2_20851.png" alt="" style="width: 49%; display: inline-block;">
</div>

![](sersicplot/24975.png)
<div>
<img src="../multicomp/img_mc_multi/24975.png" alt="" style="width: 49%; display: inline-block;">
<img src="../multicomp/img_mc_single/nircam2_24975.png" alt="" style="width: 49%; display: inline-block;">
</div>

![](sersicplot/23395.png)
<div>
<img src="../multicomp/img_mc_multi/23395.png" alt="" style="width: 49%; display: inline-block;">
<img src="../multicomp/img_mc_single/nircam3_23395.png" alt="" style="width: 49%; display: inline-block;">
</div>


### 3. Mean Mag > 26.5 mag

![](sersicplot/38262.png)
<div>
<img src="galfitmimage_linear/nircam1_38262.png" alt="" style="width: 49%; display: inline-block;">
<img src="singlebandimage/nircam1_38262.png" alt="" style="width: 49%; display: inline-block;">
</div>

![](sersicplot/41239.png)
<div>
<img src="galfitmimage_linear/nircam1_41239.png" alt="" style="width: 49%; display: inline-block;">
<img src="singlebandimage/nircam1_41239.png" alt="" style="width: 49%; display: inline-block;">
</div>
	
![](sersicplot/35552.png)
<div>
<img src="../multicomp/img_mc_multi/35552.png" alt="" style="width: 49%; display: inline-block;">
<img src="../multicomp/img_mc_single/nircam6_35552.png" alt="" style="width: 49%; display: inline-block;">
</div>


### Bad Fit or Strange Results

![](sersicplot/21844.png)
<div>
<img src="../multicomp/img_mc_multi/21844.png" alt="" style="width: 49%; display: inline-block;">
<img src="../multicomp/img_mc_single/nircam1_21844.png" alt="" style="width: 49%; display: inline-block;">
</div>
Inverse Slope

![](sersicplot/29960.png)
<div>
<img src="galfitmimage_linear/nircam1_29960.png" alt="" style="width: 49%; display: inline-block;">
<img src="singlebandimage/nircam1_29960.png" alt="" style="width: 49%; display: inline-block;">
</div>
An monotonically varying size might lead to a bad fit here. 


![](sersicplot/36301.png)
<div>
<img src="../multicomp/img_mc_multi/36301.png" alt="" style="width: 49%; display: inline-block;">
<img src="../multicomp/img_mc_single/nircam6_36301.png" alt="" style="width: 49%; display: inline-block;">
</div>


![](sersicplot/13297.png)
<div>
<img src="../multicomp/img_mc_multi/13297.png" alt="" style="width: 49%; display: inline-block;">
<img src="../multicomp/img_mc_single/nircam6_13297.png" alt="" style="width: 49%; display: inline-block;">
</div>
Sudden Jump of Sérsic Index