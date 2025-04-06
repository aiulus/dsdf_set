function addTargetLine(target, len)
    plot(0:len - 1, target * ones(1, len), 'g--', 'DisplayName', 'Target');
end