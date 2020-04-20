function [x_fit,y_fit,a,p]=level_set(x,y,G,pixel)

q=0;
p=0;%perimeter
a=0;%area

for i=2:length(y)-1
    for j=2:length(x)-1
        
        if G(i,j)==0
            q=q+1;
            m(q)=x(j);  %m =j is col
            n(q)=y(i);  %n=i is row 
        end
       if G(i,j)<0 & G(i,j+1)>0 
              q=q+1;
              m(q)=x(j)+(0-G(i,j))/(G(i,j+1)-G(i,j))*pixel;
              n(q)=y(i);
        end
        if G(i,j)>0 & G(i,j+1)<0 
              q=q+1;
              %m(q)=x(j)-gx(i,j);
              m(q)=x(j)+(G(i,j)-0)/(G(i,j)-G(i,j+1))*pixel;
              n(q)=y(i);
        end
         if G(i,j)<0 & G(i+1,j)>0 
              q=q+1;
              m(q)=x(j);
              %n(q)=y(i)+gy(i,j); 
              n(q)=y(i)+(0-G(i,j))/(G(i+1,j)-G(i,j))*pixel;
             end
          if G(i,j)>0 & G(i+1,j)<0 
              q=q+1;
              m(q)=x(j);
              %n(q)=y(i)-gy(i,j); 
              n(q)=y(i)+(G(i,j)-0)/(G(i,j)-G(i+1,j))*pixel;
              end
          if G(i,j)<0 & G(i+1,j+1)>0  
              q=q+1;
              m(q)=x(j)+(0-G(i,j))/(G(i+1,j+1)-G(i,j))*pixel;
              n(q)=y(i)+(0-G(i,j))/(G(i+1,j+1)-G(i,j))*pixel;
          end
          if G(i,j)>0 & G(i+1,j+1)<0  
              q=q+1;
              m(q)=x(j)+(G(i,j)-0)/(G(i,j)-G(i+1,j+1))*pixel;
              n(q)=y(i)+(G(i,j)-0)/(G(i,j)-G(i+1,j+1))*pixel;
          end
          if G(i,j)<0 & G(i+1,j-1)>0 
              q=q+1;
              m(q)=x(j)-(0-G(i,j))/(G(i+1,j-1)-G(i,j))*pixel;
              n(q)=y(i)+(0-G(i,j))/(G(i+1,j-1)-G(i,j))*pixel;
          end
          if G(i,j)>0 & G(i+1,j-1)<0 
              q=q+1;
              m(q)=x(j)-(G(i,j)-0)/(G(i,j)-G(i+1,j-1))*pixel;
              n(q)=y(i)+(G(i,j)-0)/(G(i,j)-G(i+1,j-1))*pixel;
          end
    end
end
x_row=m; 
y_row=n; %guoyi ser up x =col

[theta,rho] = cart2pol(x_row-mean(x_row),y_row-mean(y_row));
fitobject = fit(theta',rho','linearinterp');
% figure()
% plot(fitobject,theta',rho')
theta_fit=-pi:0.01:pi;
rho_fit = feval(fitobject,theta_fit');
theta_p=theta_fit(2)-theta_fit(1);
for i=1:length(theta_fit)
    l=rho_fit(i)*theta_p;
    p=p+l;
    a=a+l*rho_fit(i)/2;
end

[x_fit,y_fit]=pol2cart(theta_fit',rho_fit);
x_fit=x_fit+mean(x_row);
y_fit=y_fit+mean(y_row);
% figure()
% plot(x_fit,y_fit)

