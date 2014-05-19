function yp=ph_ode_generator(Lambda,theta,dim,y)
    yp=zeros(dim*(dim+2),1);
    yp(1:dim)=y(1:dim)'*Lambda;
    yp((dim+1):(2*dim))=Lambda*y((dim+1):(2*dim));
    for i=1:dim
        yp((2*dim+1+(i-1)*dim):(2*dim+i*dim))=Lambda*y((2*dim+1+(i-1)*dim):(2*dim+i*dim))+y(i)*theta;
    end
    
    
    