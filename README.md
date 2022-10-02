# Pearson & Manders
## SpatialCorrelation

![image](https://user-images.githubusercontent.com/113156852/193446180-eaa49889-225d-4344-a9b7-9f066f1880b4.png)

The digital image data were vectorized and the pixel-wise fluorograms of the markers were plotted. Then, the statistical covariance, overlap, and correlation coefficients (Pearson and Menders coefficients Equations 1-2) were calculated for each image according to the following equations. The analysis pipeline was developed in Python 3.8.
Let x be representative for A signal and y for B signal. Then we have:

$$Eq.1 \quad\quad\quad\quad\ Pearsons(x,y)= cov(x,y)/\sqrt{var(x)var(y)} $$

$$Eq.2 \quad\quad\quad\quad\ Manders(x,y)= dot(x,y)/\sqrt{dot(x,x)dot(y,y)} $$

Where cov is the covariance, var is variance, and dot represents the inner product of the variables.
The response range is between 0 and 1 for Menders and -1 to 1 for Pearson’s coefficients. The coefficient value closer to 0 shows uncorrelated data in both cases, while the coefficient closer to 1 indicates a perfect positive correlation and -1 indicates a complete negative correlation (for Pearson).  
