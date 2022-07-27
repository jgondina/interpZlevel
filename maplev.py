import numpy as np

def maplev(a):
   # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   # % This file fills the land points or invalid data points using 5-point
   # % laplacian filter. For water points(valid points), their original
   # % values are unchanged
   # %
   # % im - row of the data
   # % jm - column of the data
   # % a - data needed to process (land points are set to be NaN)
   # % out- results
   # %
   # % R.He original on 12/30/99
   # % R.He revised  on 05/06/01
   # % R.He revised  on 10/12/03
   # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   im, jm = a.shape

   ii = np.where(~np.isnan(a))
   av = a[ii].mean() # Mean value of the valid data points

   jj = np.where(np.isnan(a))
   if len(ii) == 0:
       return

   if len(ii) == 1:
       a[jj]=av
   else:
        [X,Y] = np.meshgrid(np.arange(1, jm), np.arange(1, im))
        a[jj] = np.griddata(X[ii],Y[ii],a[ii],X[jj],Y[jj],'nearest')

   b=a.copy()

   # Five point laplacian filter smoothing.
   lpp=100   # do 100 times smoothing, you may need another number. More loop, a smoother field.

   for k in range(lpp):
       i=[1:-1]
       j=[1:-1]
       cc(i,j)=b(i,j)+0.5/4*(b(i+1,j)+b(i,j-1)+b(i-1,j)+b(i,j+1)-4*b(i,j));

       # Set the boundary equal ot the next interior points
       cc(1,:) =cc(2,:);
       cc(im,:)=cc(im-1,:);
       cc(:,1) =cc(:,2);
       cc(:,jm)=cc(:,jm-1);

       b(jj)=cc(jj);

    a(jj)=cc(jj);    %% only change the invalid data points, keep valid points
                     %% unchange

   out=a;
