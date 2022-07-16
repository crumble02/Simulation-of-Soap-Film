## SOAP FILM SIMULATION (Roll: bs2016, B.stat(Hons) 1st year)

# Information about the problem:

# given boundaries: 0, 1, x, x^3
# we know, the soap film adjusts such that it attains the least possible potential energy
# if the surface is given by u(x,y), then 
# potential energy E(u) = c'Mc, where c_j = u(x,y) at the j-th vertex

# we have seen that minimising E(u) is equivalent to solving:
# M11 c1  + M12 c2 = 0, => M11 c1 = -M12 c2

# M12 needs to be multiplied to c2, which has dimension 200
# thus, M12 has 200 cols
# both have the same number of rows that equals the number of rows in c1 49^2
# therefore, M11 = 2401 * 2401, M12 = 2401 * 200

install.packages("rgl"); # installs the package
library(rgl); # we import this library to produce 3D plots

# we have partitioned c into 2x1 matrix with elements c1 and c2
# c2 contains the known frame heights (at boundary points)
# c1 contains the unknown heights (at interior points)
c2 = rep(0, 200);
x = rep(0, 200);
y = rep(0, 200);

for (i in 51 : 100){
  
  c2[i] = (i-50)/50;
  x[i] = (i-50)/50;
  y[i] = 0;
}
for (i in 101:150)
{
  c2[i] = 1;
  x[i] = 1;
  y[i] = (i-100)/50;
}
for (i in 151:200)
{
  c2[i] = ((200 - i + 1)/50)^3; # for boundary x^3
  x[i] = (200 - i + 1)/50;
  y[i] = 1;
}
for(i in 1:50)
{
  y[i] = i/50;
}

# first we find the M12 matrix
# we observe that it does not have any diagonal element
# the i,j-th element of the M12 matrix is the i , (j + 2401) th element of M
# area of each triangle is 0.0002

# the elements of M takes the values 0, -1, 4
# we store them in m
m = 0;
v = rep(0, 2401); # v will store the result at the end of each condition

for (j in 1 : 200){
  
  k <- 2401 + j;
  
  if(k == 2453){
    
    i = 49;
    m = -1;
    v[i] = v[i] + m*c2[j];
  } 
  else if(k == 2402){
    
    # only one possibilty: i = 1
    
    m = 0;
    
    v[1] = v[1] + c2[j]*m;
  } 
  else if(k < 2452){
    
    # here j is on the first row
    
    i <- k - 2402;
    m = -1;
    
    v[i] = v[i] + c2[j]*m;
    
    i <- i + 1;
    m = 0;
    
    v[i] = v[i] + c2[j]*m;
  } 
  else if(k > 2452 && k < 2502){

    if(k != 2453){
      
      i = (k-2452)*49;
      m = -1;
      
      v[i] = v[i] + c2[j]*m; 
      
      i <- i + 49;
      m = 0;
      
      v[i] = v[i] + c2[j]*m; 
    }
  } 
  else if(k >= 2502 && k < 2552){
    
    if(k == 2551){
      
      i = 2353;
      m = -1;
      
      v[i] = v[i] + c2[j]*m; 
      
    }
    else if (k == 2502){
      
      i = 2401;
      m = 0;
      
      v[i] = v[i] + c2[j]*m;
    } 
    else{
      
      i = 2401 - k + 2503;
      m = -1;
      
      v[i] = v[i] + c2[j]*m;
      
      i <- i - 1;
      m = 0;
      
      v[i] = v[i] + c2[j]*m; 
    }
  } 
  else if(k >= 2553 && k < 2601){
    
    if(k == 2553){
      
      i = 2353;
      m = -1;
      
      v[i] = v[i] + c2[j]*m; 
    }
    else if(k == 2601){
      
      i = 1;
      m = -1;
      
      v[i] = v[i] + c2[j]*m;
      
      i = 50;
      m = 0;
      
      v[i] = v[i] + c2[j]*m; 
    } 
    else{
      i = -(49*(k-2552)) + 2353; 
      m = -1;
      
      v[i] = v[i] + c2[j]*m;
      
      i = i + 49;
      m = 0;
      
      v[i] = v[i] + c2[j]*m; 
    }
  }
  
}
# here, v = M12 c2
# our objective reduces to solving: M11 c1 = -v 

# we use Gauss-Jacobi method to solve it
# we don't have to store M11
# it saves memory as M11 is a 2401x2401 sub-matrix

v = - v;
c1 = rep(0, 2401); 
# this will store the initial values for Gauss-Jacobi iteration
# we start with the guess of all being 0
x = rep(0, 2401);
y = rep(0, 2401);

for(l in 1:1000){ # we run the iterations 1000 times
  
  # this loop accounts for the Gauss-Jacobi iterations
  for(i in 1:2401){
    
    if(i == 1){
      
      j = 1; 
      m1 = 4;
      
      j = 2;
      m2 = -1;
      
      j = 50;
      m3 = -1;
      
      j = 51;
      m4 = 0;
      
      c1[1] = (v[1] - c1[2]*m2 - c1[50]*m3 - c1[51]*m4)/m1;
    } 
    else if(i>1 && i<49){
      
      j = i;
      m1 = 4;
      
      j = i - 1;
      m2 = -1;
      
      j = i + 1;
      m3 = -1;
      
      j = i + 49;
      m4 = -1;
      
      j = i + 50;
      m5 = 0;
      
      c1[i] = (v[i] - c1[i-1]*m2  - c1[i+1]*m3  - c1[i+49]*m4 - c1[i+50]*m5)/m1;
      
    } 
    else if (i == 49){
      
      j = i;
      m1 = 4;
      
      j = i - 1;
      m2 = -1;
      
      j = i * 2;
      m3 = -1;
      
      c1[i] = (v[i] - c1[i-1]*m2 - c1[i*2]*m3)/m1;
      
    } 
    else if( (i %% 49) == 0 && (i / 49) >= 2 && (i != 2401)){
      
      j = i;
      m1 = 4;
      
      j = i - 49;
      m2 = -1;
      
      j = i + 49;
      m3 = -1;
      
      j = i - 1;
      m4 = -1;
      
      j = i - 1 - 49;
      m5 = 0 ;
      
      c1[i] = (v[i] - c1[i-49]*m2 - c1[i + 49]*m3 - c1[i-1]*m4 - c1[i-50]*m5)/m1;
      
    } 
    else if(i == 2401){
      
      j = i;
      m1 = 4;
      
      j = i - 49;
      m2 = -1;
      
      j = i - 1;
      m3 = -1;
      
      j = i - 1 - 49;
      m4 = 0;
      
      c1[i] = (v[i] - c1[i-49]*m2 - c1[i-1]*m3 - c1[i-50]*m4)/m1;
    } 
    else if(i<2401 && i>2353){
      
      j = i;
      m1 = 4;
      
      j = i + 1;
      m2 = -1;
      
      j = i - 1;
      m3 = -1;
      
      j = i - 49;
      m4 = -1;
      
      j = i - 50;
      m5 = 0;
      
      c1[i] = (v[i] - c1[i+1]*m2 - c1[i-1]*m3 - c1[i-49]*m4 - c1[i-50]*m5)/m1;
    } 
    else if(i == 2353){
      j = i;
      m1 = 4;
      j = i + 1;
      m2 = -1;
      j = i - 49;
      m3 = -1;
      c1[i] = (v[i] - m2 * c1[i + 1] - m3 * c1[i - 49]) / m1;
      
    } else if( (i - 1) %% 49 == 0 && i > 1 && i < 2353)
    {
      j = i;
      m1 = 4;
      
      j = i + 1;
      m2 = -1;
      
      j = i + 49;
      m3 = -1;
      
      j = i - 49;
      m4 = -1;
      
      j = i + 1 + 49;
      m5 = 0;
      
      c1[i] = (v[i] - c1[i+1]*m2 - c1[i+49]*m3 - c1[i-49]*m4 - c1[i+50]*m5)/m1;
    } 
    else{
      
      # interior points
      j <- i;
      m1 = 4;
      
      j <- i-1; 
      m2 = -1;
      
      j <- i + 1; 
      m3 = -1;
      
      j <- i - 49; 
      m4 = -1;
      
      j <- i + 49; 
      m5 = -1;
      
      j <- i - 50; 
      m6 = 0;
      
      j <- i + 50; 
      m7 = 0;
      
      c1[i] = (v[i] - c1[i-1]*m2 - c1[i+1]*m3 - c1[i-49]*m4 - c1[i+49]*m5 - c1[i-50]*m6 - c1[i+50]*m7)/m1;
    }
    x[i] = (49 - i %% 49)/50;
    y[i] = (49 - (i/49))/50;
    
  }
}

plot3d(x, y, c1, size = 3, col="purple", axes = TRUE);
# this line plots the soap film in 3D