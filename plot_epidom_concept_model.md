```python
import numpy as np
import matplotlib as mpl
import os
import pandas as pd
```


```python
import matplotlib.pyplot as plt
```


```python
import seaborn as sns
```


```python
sns.set_style("whitegrid")
sns.color_palette("colorblind")
sns.set(font_scale=1)

plt.rcParams["font.family"] = "Lucida Grande"
plt.rcParams["font.size"] = 12

```


```python
#make this example reproducible
np.random.seed(1)

#generate array of 1000 values that follow uniform distribution with min=0 and max=1

######FULL ADDITIVITY, h= 0.5, allele substitution effect = a*0.5 (as d = 0, h = 0.5) #######
##allele substitution effect alpha = ah - 2dq

p = np.random.uniform(low=0, high=1, size=1000)
q = 1-p
a = 10
h = 0.5
d = a*(h-0.5)
alpha = a*h-d*2*(1-p)
```


```python
sigmaA = 2*p*q*alpha*alpha
sigmaD = 4*d*d*p*p*q*q
```


```python
data = {'p' : p, 'sigmaA' : sigmaA, 'sigmaD' : sigmaD, 'alpha' : alpha}
df = pd.DataFrame(data)
```


```python
sns.set_style("whitegrid")
sns.color_palette("colorblind")
sns.set(font_scale=1)

plt.rcParams["font.family"] = "Lucida Grande"
plt.rcParams["font.size"] = 12
```


```python
sns.set_style("whitegrid")

fig1 = sns.lmplot(x='p', y='value', hue = 'variable',
           data = pd.melt(df, ['p']), height =5, aspect = 2,
           order = 2, scatter_kws={"s": 5}).set(title = "No dominance h = 0.5")
```


    
![png](output_8_0.png)
    



```python
##FULL dominance, h= 1, 
#allele substitution effect alpha =  (ah-2dq)
```


```python
p = np.random.uniform(low=0, high=1, size=1000)
q = 1-p
a = 10
h = 1.0
d = a*(h-0.5)
alpha = a*h-d*2*(1-p)
```


```python
sigmaA = 2*p*q*alpha*alpha
sigmaD = 4*d*d*p*p*q*q
```


```python
data = {'p' : p, 'sigmaA' : sigmaA, 'sigmaD' : sigmaD, 'alpha' : alpha}
df = pd.DataFrame(data)
```


```python
sns.set_style("whitegrid")

fig2 = sns.lmplot(x='p', y='value', hue = 'variable',
           data = pd.melt(df, ['p']), height =5, aspect = 2,
           order = 4, scatter_kws={"s": 5}).set(title = "Full dominance h = 1")
```


    
![png](output_13_0.png)
    



```python
###IF h = 0.75

p = np.random.uniform(low=0, high=1, size=1000)
q = 1-p
a = 10
h = 0.75
d = a*(h-0.5)
alpha = a*h-d*2*(1-p)

sigmaA = 2*p*q*alpha*alpha
sigmaD = 4*d*d*p*p*q*q


data = {'p' : p, 'sigmaA' : sigmaA, 'sigmaD' : sigmaD, 'alpha' : alpha}
df = pd.DataFrame(data)


sns.set_style("whitegrid")

fig3 = sns.lmplot(x='p', y='value', hue = 'variable',
           data = pd.melt(df, ['p']), height =5, aspect = 2,
           order = 4, scatter_kws={"s": 5}).set(title = "h = 0.75")
```


    
![png](output_14_0.png)
    



```python
###IF h = 0.25

p = np.random.uniform(low=0, high=1, size=1000)
q = 1-p
a = 10
h = 0.25
d = a*(h-0.5)
alpha = a*h-d*2*(1-p)

sigmaA = 2*p*q*alpha*alpha
sigmaD = 4*d*d*p*p*q*q


data = {'p' : p, 'sigmaA' : sigmaA, 'sigmaD' : sigmaD, 'alpha' : alpha}
df = pd.DataFrame(data)


sns.set_style("whitegrid")
fig4 = sns.lmplot(x='p', y='value', hue = 'variable',
           data = pd.melt(df, ['p']), height =5, aspect = 2,
           order = 4, scatter_kws={"s": 5}, palette = "colorblind").set(title = "h = 0.25")
```


    
![png](output_15_0.png)
    



```python
###IF h = 0

p = np.random.uniform(low=0, high=1, size=1000)
q = 1-p
a = 10
h = 0
d = a*(h-0.5)
alpha = a*h-d*2*(1-p)

sigmaA = 2*p*q*alpha*alpha
sigmaD = 4*d*d*p*p*q*q


data = {'p' : p, 'sigmaA' : sigmaA, 'sigmaD' : sigmaD, 'alpha' : alpha}
df = pd.DataFrame(data)


sns.set_style("whitegrid")
fig5 = sns.lmplot(x='p', y='value', hue = 'variable',
           data = pd.melt(df, ['p']), height =5, aspect = 2,
           order = 4, scatter_kws={"s": 5}, palette = "colorblind").set(title = "h = 0")
```


    
![png](output_16_0.png)
    



```python
#Under-dominance h < -1

p = np.random.uniform(low=0, high=1, size=1000)
q = 1-p
a = 10
h = -2
d = a*(h-0.5)
alpha = a*h-d*2*(1-p)

sigmaA = 2*p*q*alpha*alpha
sigmaD = 4*d*d*p*p*q*q


data = {'p' : p, 'sigmaA' : sigmaA, 'sigmaD' : sigmaD, 'alpha' : alpha}
df = pd.DataFrame(data)


sns.set_style("whitegrid")
fig6 = sns.lmplot(x='p', y='value', hue = 'variable',
           data = pd.melt(df, ['p']), height =5, aspect = 2,
           order = 4, scatter_kws={"s": 5}, palette = "colorblind").set(title = "Underdominance h = -2")

```


    
![png](output_17_0.png)
    



```python
#Over dominance h > 1

p = np.random.uniform(low=0, high=1, size=1000)
q = 1-p
a = 10
h = 2
d = a*(h-0.5)
alpha = a*h-d*2*(1-p)

sigmaA = 2*p*q*alpha*alpha
sigmaD = 4*d*d*p*p*q*q


data = {'p' : p, 'sigmaA' : sigmaA, 'sigmaD' : sigmaD, 'alpha' : alpha}
df = pd.DataFrame(data)


sns.set_style("whitegrid")
fig7 = sns.lmplot(x='p', y='value', hue = 'variable',
           data = pd.melt(df, ['p']), height =5, aspect = 2,
           order = 4, scatter_kws={"s": 5}, palette = "colorblind").set(title = "Overdominance h = 2")

```


    
![png](output_18_0.png)
    



```python
os.chdir('/Users/nhutran/Documents/PhD/epidom/model_plots/Epistasis')
```


```python
fig7.savefig('Additive_var_over_dominance.pdf')
```


```python
###############################################################
```


```python
###############################################################
```


```python
#####################EPISTATIC MODEL 1 #####################
```


```python
#for 2 loci A and B
#A has p(A1) and q(A2) allele frequencies
#B has r(B1) and s(B2) allele frequencies
#epistastic factor k 
    # when k > 1, we speak of synergistic epistasis: sum of double homzygotes is greater than separate ones
    # when k < 1, we speak of antagonistic epistasis
    # when k = 0, there is no epistasis
     #f.e.here, value in excess of A2A2B2B2 is c = (k-1)(a(A) + a(B))
# call 'aA' and 'aB' as additive value of A2A2 and B2B2 respectively
# call 'hA' and 'hB' as dominance coefficient of locus A and B respectively


#differnt dominance coefficients
 #1. as textbook: hA = 0.8, hB = 0.5 # so this is somewhat dominance at one locus and no dominance at the other
 #2. no dominance at both locus h = 0.5
 #3. full dominance at both locus h = 1
 #4. full dominance at one locus, no dominance at the other # pretty much like text book
 #5. somewhat dominance at both locus,say 0.7 and 0.8


#locus A
aA = 10
hA = 0.7
p = np.random.uniform(low=0, high=1, size=1000)
q = 1 - p 

#locus B
aB = 12
hB = 0.8
r = 0.1
s = 1 - r

#excess value of double homozygotes due to epistasis
k = 2
c = (k-1)*(aA + aB)

#so we have at locus A
dA = aA*(hA-0.5)
dAe = dA - (0.5*c)*s*s #dominance value due to epistasis
alphaAe = aA*hA - 2*dAe*q #allele substitution due to epistasis

#and at locus B
dB = aB*(hB-0.5)
dBe = dB - (0.5*c)*q*q
alphaBe = aB*hB - 2*dBe*s


#epistatic effects

alphaAe_alphaBe = c*q*s
alphaAe_dBe = (0.5*c)*q
dAe_alphaBe = (0.5*c)*s
dAe_dBe = 0.25*c

#Then global components will be

M = (aA*q + 2*dA*p*q) + (aB*s + 2*dB*r*s) + c*q*q*s*s # population mean
Va = 2*(alphaAe*alphaAe)*p*q + 2*(alphaBe*alphaBe)*r*s
Vd = ((2*dAe*p*q)**2) + ((2*dBe*r*s)**2)
Vaa = 4*((c*q*s)**2)*p*q*r*s
Vad = 2*((c*q*s)**2)*(p*q*r*r + p*p*r*s)
Vdd = (c*p*q*r*s)**2
Vepi = Vaa + Vad + Vdd 
Vg = Va + Vd + Vepi

```


```python
data = {'p' : p, 'Va' : Va, 'Vd' : Vd, 'Vepi' : Vepi, 'alphaA' : alphaAe, 'alphaB': alphaBe}
df = pd.DataFrame(data)

```


```python
sns.set_style("whitegrid")
fig8 = sns.lmplot(x='p', y='value', hue = 'variable',
           data = pd.melt(df, ['p']), height =5, aspect = 2,
           order = 4, scatter_kws={"s": 3}, 
                  palette = "colorblind").set(title = "hA = 0.7, hB = 0.8, r = 0.1")

```


    
![png](output_26_0.png)
    



```python
sns.set_style("whitegrid")
fig12 = sns.lmplot(x='p', y='value', hue = 'variable',
           data = pd.melt(df, ['p']), height =5, aspect = 2,
           order = 4, scatter_kws={"s": 3}, 
                   palette = "colorblind").set(title = "hA = 0.7, hB = 0.8, r = 0.5")

```


    
![png](output_27_0.png)
    



```python
sns.set_style("whitegrid")
fig10 = sns.lmplot(x='p', y='value', hue = 'variable',
           data = pd.melt(df, ['p']), height =5, aspect = 2,
           order = 4, scatter_kws={"s": 3}, 
                   palette = "colorblind").set(title = "hA = 0.7, hB = 0.8, r = 0.9")

```


    
![png](output_28_0.png)
    



```python
fig8.savefig('Vepi_high_recessive_hA0.7_hB0.8.pdf')
```


```python
fig12.savefig('Vepi_mid_freq_hA0.7_hB0.8.pdf')
```


```python
fig10.savefig('Vepi_low_recessive_hA0.7_hB0.8.pdf')
```


```python
#########################################################################################################
```


```python
#################### EPISTATIC MODEL 2 - trying to make bar plots like in the book #####################
#####################             but did not look so nice       #######################################
```


```python
#for 2 loci A and B
#A has p(A1) and q(A2) allele frequencies
#B has r(B1) and s(B2) allele frequencies
#epistastic factor k 
    # when k > 1, we speak of synergistic epistasis: sum of double homzygotes is greater than separate ones
    # when k < 1, we speak of antagonistic epistasis
    # when k = 0, there is no epistasis
     #f.e.here, value in excess of A2A2B2B2 is c = (k-1)(a(A) + a(B))
# call 'aA' and 'aB' as additive value of A2A2 and B2B2 respectively
# call 'hA' and 'hB' as dominance coefficient of locus A and B respectively


#locus A
aA = 10
hA = 0.8
p = np.array([0.1, 0.5, 0.9])
q = 1 - p 

#locus B
aB = 12
hB = 0.5
r = 0.9
s = 1 - r

#excess value of double homozygotes due to epistasis
k = 2
c = (k-1)*(aA + aB)

#so we have at locus A
dA = a*(h-0.5)
dAe = dA - (0.5*c)*s*s #dominance value due to epistasis
alphaAe = aA*h - 2*dAe*q #additive values due to epistasis

#and at locus B
dB = aB*(hB-0.5)
dBe = dB - (0.5*c)*q*q
alphaBe = aB*hB - 2*dBe*s


#epistatic effects

alphaAe_alphaBe = c*q*s
alphaAe_dBe = (0.5*c)*q
dAe_alphaBe = (0.5*c)*s
dAe_dBe = 0.25*c

#Then global components will be

M = (aA*q + 2*dA*p*q) + (aB*s + 2*dB*r*s) + c*q*q*s*s # population mean
Va = 2*(alphaAe*alphaAe)*p*q + 2*(alphaBe*alphaBe)*r*s
Vd = ((2*dAe*p*q)**2) + ((2*dBe*r*s)**2)
Vaa = 4*((c*q*s)**2)*p*q*r*s
Vad = 2*((c*q*s)**2)*(p*q*r*r + p*p*r*s)
Vdd = (c*p*q*r*s)**2
Vepi = Vaa + Vad + Vdd 
Vg = Va + Vd + Vepi

```


```python
Additive = Va*100/Vg
Dominance = Vd*100/Vg
Epistasis = Vepi*100/Vg
```


```python
r1 = pd.DataFrame({'p' : p, 'r' : r,'Va' : Additive, 'Vd' : Dominance, 'Vi' : Epistasis})

```


```python
r1.dtypes
```




    p     float64
    r     float64
    Va    float64
    Vd    float64
    Vi    float64
    dtype: object




```python
r1['r'] = r1['r'].astype('category')
r1['p'] = r1['p'].astype('category')
```


```python
r11 = pd.melt(r1, id_vars = ['p'], value_vars=['Va','Vd','Vi'])
```


```python
r11.dtypes
```




    p           category
    variable      object
    value        float64
    dtype: object




```python
sns.set_style("whitegrid")
sns.catplot(data=r11, x="p", y="value", hue="variable", kind="bar").set(title = "r = 0.9") # r = 0.9
```

    /Users/nhutran/opt/anaconda3/lib/python3.9/site-packages/pandas/io/formats/format.py:1429: FutureWarning:
    
    Index.ravel returning ndarray is deprecated; in a future version this will return a view on self.
    





    <seaborn.axisgrid.FacetGrid at 0x7fd21a5319d0>




    
![png](output_42_2.png)
    



```python
sns.set_style("whitegrid")
sns.catplot(data=r11, x="p", y="value", hue="variable", kind="bar").set(title = "r = 0.5") # r = 0.5
```

    /Users/nhutran/opt/anaconda3/lib/python3.9/site-packages/pandas/io/formats/format.py:1429: FutureWarning:
    
    Index.ravel returning ndarray is deprecated; in a future version this will return a view on self.
    





    <seaborn.axisgrid.FacetGrid at 0x7fd23def7970>




    
![png](output_43_2.png)
    



```python
sns.set_style("whitegrid")
sns.catplot(data=r11, x="p", y="value", hue="variable", kind="bar").set(title = "r = 0.1") # r = 0.1
```

    /Users/nhutran/opt/anaconda3/lib/python3.9/site-packages/pandas/io/formats/format.py:1429: FutureWarning:
    
    Index.ravel returning ndarray is deprecated; in a future version this will return a view on self.
    





    <seaborn.axisgrid.FacetGrid at 0x7fd1eb7448b0>




    
![png](output_44_2.png)
    

