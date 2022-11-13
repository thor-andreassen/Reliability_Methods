function [Zamv]=AMV(func_handle,models,mus,P,tol,use_approx)   
    X=mus;
    for count_var=1:length(X)
            [X_means_RF(count_var),X_stds_RF(count_var)]=rackwitzFiesslerApprox(models(count_var),X(count_var));
    end
    mus=X_means_RF';
    sigmas=X_stds_RF';
    a0=func_handle(mus);
    [dz_dx]=numericGradient(func_handle, mus);
    
    dz_dx_norm=dz_dx.*sigmas;
        
    %% Determine the Alpha/Direction Cosines
        alphas=dz_dx_norm/norm(dz_dx_norm);
        
        %% Determine the Output Value
        Zamv.Z=zeros(size(P));
        for i=1:length(P)
            beta=norminv(P(i));
            Xmpp=mus+(beta.*alphas).*sigmas;
            if use_approx
                Zamv.Z(i)=a0+sum(dz_dx.*Xmpp);
            else
                Zamv.Z(i)=func_handle(Xmpp);
            end
            
        end
        Zamv.probabilities=P;
    end