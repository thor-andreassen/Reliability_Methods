function plotTrendLine(x,y,order)
    p = polyfit(x, y, order);
    px = [min(x) max(x)];
    py = polyval(p, px);
    scatter(x, y,'.')
    hold on
    plot(px, py, 'LineWidth', 2);
    str_val='';
    flip_p=fliplr(p);
    for count_order=1:length(p)
        str_val=[str_val,num2str(flip_p(count_order)),'x^',num2str(count_order-1)];
        if count_order<length(p)
            str_val=[str_val,' + '];
        end
    end
    title(str_val)
end