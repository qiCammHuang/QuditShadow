function Variance = get_stat_variance(data)

val_mean = mean(data);

Variance = 0;
for i = 1:length(data)
    Variance = Variance + (data(i) - val_mean)^2;
end

Variance = sqrt(Variance/(length(data)));

end