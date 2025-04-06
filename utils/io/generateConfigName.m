function configname = generateConfigName(descriptions, runIdx)
    prefix = string(descriptions(runIdx));
    suffix = sprintf('-run%d', runIdx);
    configname = strcat(prefix, suffix);
end