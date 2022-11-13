function [norm_mean_approx,norm_std_approx]=rackwitzFiesslerApprox(cdf_model,val)
    norm_std_approx=zeros(length(val),1);
    norm_mean_approx=zeros(length(val),1);
    for counti=1:length(val)
        Fx_val=cdf(cdf_model,val(counti));
        fx_val=pdf(cdf_model,val(counti));
        norm_std_approx(counti)=(normpdf(norminv(Fx_val)))/(fx_val);
        norm_mean_approx(counti)=val(counti)-norm_std_approx(counti)*norminv(Fx_val);
    end
end