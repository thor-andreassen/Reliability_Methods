function [Beta,Pf,alphas,X]=form2(func_handle,models,mus,limit_state,tol)
    X_0=mus;
    limit_func_handle=@(in1) limit_state-func_handle(in1);
    beta_tol=tol;
    g_tol=tol;
    
    conv_g=Inf;
    conv_beta=Inf;
    max_count=10;
    
    Beta_old=Inf;
    
    counter=1;
    
    X=X_0;
    X_means_RF=zeros(length(X),1);
    X_stds_RF=zeros(length(X),1);
    while counter<max_count && (conv_beta>beta_tol || conv_g>g_tol)
        % RF approximation of mean and standard deviation for non normal
        % variables
        for count_var=1:length(X)
            [X_means_RF(count_var),X_stds_RF(count_var)]=rackwitzFiesslerApprox(models(count_var),X(count_var));
        end
        
        X_star=(X-X_means_RF)./X_stds_RF;
        func_val=limit_func_handle(X);
        grad_val=numericGradient(limit_func_handle,X);
        
        grad_val_star=grad_val.*X_stds_RF;
        
        alphas=grad_val_star/norm(grad_val_star);
        
        X_star=(1/norm(grad_val_star)^2)*(grad_val_star'*X_star-func_val)*grad_val_star;
        
        Beta=norm(X_star);
        conv_beta=abs(Beta-Beta_old);
        Beta_old=Beta;
        conv_g=abs(func_val);
        
        X=X_means_RF+X_star.*X_stds_RF;
        
        counter=counter+1;
    end
    Pf=1-normcdf(Beta);
end