function coeff=lf(x,y,lb,ub)

xmax=max(x);
ymax=max(y);

for i=1:size(x,2)
    xnorm(:,i)=x(:,i)/max(1e-3,xmax(i));   
    ynorm(:,i)=y(:,i)/max(1e-3,ymax(i)); 
end

val = xnorm+ynorm;

val_ub = prctile(val, ub);
val_lb = prctile(val, lb);
indi=[];

for j=1:size(val,1)
    if (val(j)>val_ub || val(j)<val_lb ) && xnorm(j)>0 && ynorm(j)>0
        indi=[indi,j];
    end
end

x_fit=x(indi,:);
y_fit=y(indi,:);

coeff=y_fit'*x_fit/norm(x_fit,2);
end