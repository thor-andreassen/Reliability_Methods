function [grad]=numericGradient(func, value, perturb)
    if nargin<3
        perturb=1E-5;
    end
    func_val=func(value);
    grad=zeros(length(value),1);
    for count_var=1:length(value)
        current_value=value;
        current_value(count_var)=current_value(count_var)+perturb;
        grad(count_var)=(func(current_value)-func_val)/perturb;
    end
    
end