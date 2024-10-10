function fjb
suh_pipelines;
if ~isdeployed 
    try
        PyEnvironment.Setup;
    catch ex
        disp(ex.message);
    end
end
end